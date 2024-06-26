﻿import os
import time

import pandas as pd
from astropy.time import Time
import utils
from sqlalchemy import text
from Astrometry import DoAss
from photometric_calibration import apply_both, make_master_dark, make_master_flat

time_wait_sec = 10
# DARK_ROOT = r'D:\_Darks\mdark_'   # mflat_exp_temp
# FLAT_ROOT = r'D:\_Flats\mflat_'
ID = 'frame_id'
PATH_NAME = 'frame_path'
CAL_PATH_NAME = 'calibration_frame_path'
DO_DARK_NAME = 'is_do_dark'
DO_FLAT_NAME = 'is_do_flat'
DO_ASTROMETRY_NAME = 'is_do_astrometry'
DO_SEX_NAME = 'is_do_sex'
TROUBLES_NAME = 'trouble_count'
# names = [PATH_NAME, DO_DARK_NAME, DO_FLAT_NAME, DO_ASTROMETRY_NAME, TROUBLES_NAME]
API_KEY = 'hipfhzhlzygnlvix'
delta_temp = 1


def sisyphus():
    print(f'\nNew sisyphus tic at {Time.now()}')
    # read bd
    eng = utils.connect_to_db()
    print('Connected to robophot db')

    query_for_current_task = ("SELECT frame_id, frame_path, calibration_frame_path, coord2000, "
                              "frame_filter, date_utc, ccd_temp, is_do_dark, is_do_flat, is_do_astrometry, "
                              "is_do_sex, trouble_count "
                              "FROM robophot_frames, robophot_tasks "
                              "WHERE fk_task_id = task_id AND (frame_type = 'Object') AND "  # OR frame_type = 'Test'
                              "status = 1 AND "
                              "(NOT is_do_dark OR NOT is_do_flat OR NOT is_do_astrometry OR NOT is_do_sex) AND "
                              "trouble_count < 3 "
                              "ORDER BY (CASE "
                              "WHEN frame_filter = 'i' THEN 1 "
                              "WHEN frame_filter = 'r' THEN 2 "
                              "WHEN frame_filter = 'g' THEN 3 END), frame_id")
    table = pd.read_sql_query(query_for_current_task, eng)
    # query = "SELECT get_frames_for_cal() FROM robophot_frames"
    if table.shape[0] == 0:
        query = ("SELECT frame_id, frame_path, calibration_frame_path, coord2000, "
                 "frame_filter, date_utc, ccd_temp, is_do_dark, is_do_flat, is_do_astrometry, is_do_sex, trouble_count "
                 "FROM robophot_frames, robophot_tasks "
                 "WHERE fk_task_id = task_id AND (frame_type = 'Object') AND "  # OR frame_type = 'Test'
                 "(NOT is_do_dark OR NOT is_do_flat OR NOT is_do_astrometry OR NOT is_do_sex) AND trouble_count < 3 "
                 "ORDER BY (CASE "
                 "WHEN frame_filter = 'i' THEN 1 "
                 "WHEN frame_filter = 'r' THEN 2 "
                 "WHEN frame_filter = 'g' THEN 3 END), frame_id")
        table = pd.read_sql_query(query, eng)

    # apply calibration
    if table.shape[0] > 0:
        print(f'Start calibrate, got {table.shape[0]} frames')
        for indx, row in table.iterrows():
            done_cal = 0
            med_fwhm = 0
            med_ell = 0
            med_bkg = 0
            dark_id = 0
            flat_id = 0
            if not row[DO_DARK_NAME]:
                q_get_for_do_dark = (f"SELECT m_frame_id, m_frame_path FROM "
                                     f"robophot_master_frames, robophot_frames, robophot_tasks "
                                     f"WHERE {row[ID]} = frame_id AND robophot_frames.fk_task_id = task_id AND "
                                     f"m_frame_type = 'm_Dark' AND m_frame_path IS NOT NULL AND "
                                     f"ABS(ccd_temp-ccd_temp_mean) < {delta_temp} AND camera_sn = m_camera_sn AND "
                                     f"x_bin = m_x_bin AND y_bin = m_y_bin AND exp_time = m_exp_time AND "
                                     f"(date_make_utc - date_utc)::INTERVAL < '2 DAY'::interval "
                                     f"ORDER BY ABS(EXTRACT(DAY FROM (date_make_utc - date_utc)::INTERVAL)) "
                                     f"LIMIT 4")
                darks_table = pd.read_sql_query(q_get_for_do_dark, eng)
                if darks_table.shape[0] == 0:
                    print(f"Can\'t calibrate frame #{row[ID]}, there are no any darks")
                    row[TROUBLES_NAME] += 1
                q_get_for_do_flat = (
                    f"SELECT m_frame_id, m_frame_path FROM robophot_master_frames, robophot_frames, robophot_tasks "
                    f"WHERE {row[ID]} = frame_id AND robophot_frames.fk_task_id = task_id AND "
                    f"m_frame_type = 'm_Flat' AND m_frame_path IS NOT NULL AND "
                    f"m_frame_filter = frame_filter AND "
                    f"camera_sn = m_camera_sn AND x_bin = m_x_bin AND y_bin = m_y_bin AND "
                    f"(date_make_utc - date_utc)::INTERVAL < '2 DAY'::interval "
                    f"ORDER BY ABS(EXTRACT(DAY FROM (date_make_utc - date_utc)::INTERVAL)) "
                    f"LIMIT 4")
                flats_table = pd.read_sql_query(q_get_for_do_flat, eng)
                if flats_table.shape[0] == 0:
                    print(f"Can\'t calibrate frame #{row[ID]} {row[PATH_NAME]}, there are no any flats")
                    row[TROUBLES_NAME] += 1
                for d_index, dark in darks_table.iterrows():
                    if row[DO_DARK_NAME] is True:
                        break
                    for f_indx, flat in flats_table.iterrows():
                        if apply_both(row[PATH_NAME], dark['m_frame_path'], flat['m_frame_path'], row[CAL_PATH_NAME]):
                            row[DO_DARK_NAME] = True
                            row[DO_FLAT_NAME] = True
                            dark_id = dark['m_frame_id']
                            flat_id = flat['m_frame_id']
                            row[TROUBLES_NAME] = 0
                            print(f'Dark and flat applied to frame #{row[ID]}, path to cal {row[CAL_PATH_NAME]}')
                            # print(f'Flat applied to frame #{row[ID]}, path {row[CAL_PATH_NAME]}')
                            break
                        else:
                            print(f"Can\'t apply dark {dark['m_frame_path']} or flat {flat['m_frame_path']} on "
                                  f"frame #{row[ID]}, path {row[CAL_PATH_NAME]}")
                        row[TROUBLES_NAME] += 1
                done_cal += 1
            if row[DO_DARK_NAME] and not row[DO_SEX_NAME]:
                med_fwhm, med_ell, med_bkg = utils.get_fwhm_data(row[CAL_PATH_NAME])
                if med_fwhm > 0:
                    print('Made sex')
                    row[DO_SEX_NAME] = True
                    row[TROUBLES_NAME] = 0
                else:
                    print(f'Can\'t do sex on frame #{row[ID]}, path {row[CAL_PATH_NAME]}')
                    row[TROUBLES_NAME] += 1
                done_cal += 1
            if row[DO_SEX_NAME] and not row[DO_ASTROMETRY_NAME]:
                if DoAss(row[CAL_PATH_NAME]):
                    row[DO_ASTROMETRY_NAME] = True
                    row[TROUBLES_NAME] = 0
                    print(f'Made WCS on frame #{row[ID]}, path {row[CAL_PATH_NAME]}')
                else:
                    row[DO_ASTROMETRY_NAME] = False
                    row[TROUBLES_NAME] += 1
                    print(f"Can\'t make WCS on frame #{row[ID]}, path {row[CAL_PATH_NAME]}")
                done_cal += 1

            # update
            if done_cal > 0:
                q_update = (f"UPDATE robophot_frames SET ({'fk_master_dark_id, ' if dark_id > 0 else ''}"
                            f"{'fk_master_flat_id, ' if dark_id > 0 else ''}"
                            f"sex_fwhm, sex_ell, sex_background, "
                            f"{DO_DARK_NAME}, {DO_FLAT_NAME}, {DO_ASTROMETRY_NAME}, {DO_SEX_NAME}, {TROUBLES_NAME}) = "
                            f"({str(dark_id) + ', ' if dark_id > 0  else ''}"
                            f"{str(flat_id) + ', ' if flat_id > 0 else ''} "
                            f"{med_fwhm}, {med_ell}, {med_bkg}, "
                            f"{row[DO_DARK_NAME]}, {row[DO_FLAT_NAME]}, {row[DO_ASTROMETRY_NAME]}, "
                            f"{row[DO_SEX_NAME]}, {row[TROUBLES_NAME]}) WHERE {ID} = {row[ID]}")
                with eng.begin() as conn:     # TRANSACTION
                    conn.execute(text(q_update))
                print()
    else:
        print('No frames to calibrate')

    # make master dark
    q_get_m_dark_frame = ("SELECT m_frame_id, m_frame_type, m_frame_filter, "
                          "m_camera_sn, m_x_bin, m_y_bin, m_exp_time, fk_task_id, m_trouble FROM "
                          "robophot_master_frames WHERE "
                          "date_make_utc IS NULL AND m_frame_type = 'm_Dark' AND m_trouble < 3")
    m_darks = pd.read_sql_query(q_get_m_dark_frame, eng)
    print(f'Start make master darks, there is {m_darks.shape[0]} darks to create')
    '''
    - У кадров должны быть одинаковыми:
      -- IMAGETYP
      -- XBINNING
      -- YBINNING
      -- SERNUM
      -- EXPTIME
    '''
    for m_d_indx, m_dark in m_darks.iterrows():
        print("D"*42)
        print(f"Make dark: working on frame #{m_dark['m_frame_id']}")
        q_get_mean_ccd_temp = ("SELECT round(avg(ccd_temp)::numeric, 1) FROM robophot_frames, robophot_tasks WHERE "
                               f"frame_type = 'Dark' AND "
                               f"fk_task_id = task_id AND task_id = '{m_dark['fk_task_id']}' AND "
                               f"{m_dark['m_x_bin']} = x_bin AND {m_dark['m_y_bin']} = y_bin AND "
                               f"'{m_dark['m_camera_sn']}' = camera_sn AND {m_dark['m_exp_time']} = exp_time")
        mean_ccd_temp = pd.read_sql_query(q_get_mean_ccd_temp, eng)['round'][0]
        set_m_trouble = ("UPDATE robophot_master_frames SET "
                         "m_trouble = 3 WHERE "
                         f"m_frame_id = {m_dark['m_frame_id']}")
        if mean_ccd_temp is None:
            print("mean_ccd_temp is None")
            with eng.begin() as conn:     # TRANSACTION
                conn.execute(text(set_m_trouble))
            continue
        q_get_darks_for_m = (f"SELECT frame_path, date_utc FROM "
                             f"robophot_frames, robophot_tasks "
                             f"WHERE frame_type = 'Dark' AND "
                             f"fk_task_id = task_id AND task_id = '{m_dark['fk_task_id']}' AND "
                             f"{m_dark['m_x_bin']} = x_bin AND {m_dark['m_y_bin']} = y_bin AND "
                             f"'{m_dark['m_camera_sn']}' = camera_sn AND {m_dark['m_exp_time']} = exp_time AND "
                             f"abs(({mean_ccd_temp} - robophot_frames.ccd_temp)::numeric) < {delta_temp}")
        '''
        - Оставляешь только те кадры, где CCD-TEMP
        отличается по модулю от mean_ccd_temp
        не более чем на epsilon_temp
        '''
        darks_for_master = pd.read_sql_query(q_get_darks_for_m, eng)
        '''
        - Если после этого в наборе осталось меньше 
        min_images_num (у меня это 5 шт),
        то создание мастер-дарка провалено
        '''
        if darks_for_master.shape[0] < 5:
            print("Too few frames to make master dark")
            with eng.begin() as conn:     # TRANSACTION
                conn.execute(text(set_m_trouble))
            continue
        date_obs = Time(darks_for_master['date_utc'][0])
        date_made = date_obs.to_value(format='fits')
        # path = DARK_ROOT + f'{m_dark["m_exp_time"]}_{mean_ccd_temp}.fits.gz'
        path_root = rf'D:\{date_obs.to_value("ymdhms")[0]}\{date_obs.to_value("iso", subfmt="date")}\UNKNOWN\MASTER_DARK'
        if not os.path.exists(path_root):
            os.makedirs(path_root)
        path = rf'{path_root}\MDARK_{date_made.replace(":", "-")}_{m_dark["m_frame_filter"]}_{m_dark["m_exp_time"]}_{mean_ccd_temp}.fits.gz'
        if make_master_dark(darks_for_master['frame_path'], mean_ccd_temp, Time.now().to_value(format='fits'), path):
            q_update_m_dark = ("UPDATE robophot_master_frames SET (m_frame_path, ccd_temp_mean, date_make_utc) = "
                               f"('{path}'::text, {mean_ccd_temp}::numeric, "
                               f"'{date_made}'::timestamp) WHERE "
                               f"m_frame_id = {m_dark['m_frame_id']}")
            with eng.begin() as conn:     # TRANSACTION
                conn.execute(text(q_update_m_dark))
            # con.execute(q_update_m_dark)
            # eng.execute(text(q_update_m_dark))
            print(f"Made master dark {path}")
        else:
            print("Failed while make master dark")
            with eng.begin() as conn:     # TRANSACTION
                conn.execute(text(set_m_trouble))

    # make flat
    q_get_m_flat_frame = ("SELECT m_frame_id, m_frame_type, m_frame_filter, "
                          "m_camera_sn, m_x_bin, m_y_bin, m_exp_time, m_trouble, fk_task_id FROM "
                          "robophot_master_frames WHERE "
                          "date_make_utc IS NULL AND m_frame_type = 'm_Flat' AND m_trouble < 3")
    m_flats = pd.read_sql_query(q_get_m_flat_frame, eng)
    print(f'Start make master flats, there is {m_flats.shape[0]} flats to create')
    if m_flats.shape[0] != 0:

        '''
        - У кадров должны быть одинаковыми:
          -- IMAGETYP
          -- XBINNING
          -- YBINNING
          -- SERNUM
          -- EXPTIME
        '''
        for m_flat_indx, m_flat in m_flats.iterrows():
            print("F"*42)
            print(f"Make flat: working on frame #{m_flat['m_frame_id']}")
            q_get_mean_ccd_temp = ("SELECT round(avg(ccd_temp)::numeric, 1) FROM robophot_frames, robophot_tasks WHERE "
                                   f"frame_type = 'Flat' AND "
                                   f"'{m_flat['m_frame_filter']}' = frame_filter AND "
                                   f"fk_task_id = task_id AND task_id = '{m_flat['fk_task_id']}' AND "
                                   f"{m_flat['m_x_bin']} = x_bin AND {m_flat['m_y_bin']} = y_bin AND "
                                   f"'{m_flat['m_camera_sn']}' = camera_sn AND {m_flat['m_exp_time']} = exp_time")
            mean_ccd_temp = pd.read_sql_query(q_get_mean_ccd_temp, eng)['round'][0]
            set_m_trouble = ("UPDATE robophot_master_frames SET "
                             "m_trouble = 3 WHERE "
                             f"m_frame_id = {m_flat['m_frame_id']}")
            if mean_ccd_temp is None:
                print("mean_ccd_temp is None")
                with eng.begin() as conn:     # TRANSACTION
                    conn.execute(text(set_m_trouble))
                continue
            q_get_flats_for_m = (f"SELECT frame_path, date_utc FROM "
                                 f"robophot_frames, robophot_tasks "
                                 f"WHERE frame_type = 'Flat' AND "
                                 f"'{m_flat['m_frame_filter']}' = frame_filter AND "
                                 f"fk_task_id = task_id AND task_id = '{m_flat['fk_task_id']}' AND "
                                 f"{m_flat['m_x_bin']} = x_bin AND {m_flat['m_y_bin']} = y_bin AND "
                                 f"'{m_flat['m_camera_sn']}' = camera_sn AND {m_flat['m_exp_time']} = exp_time AND "
                                 f"abs(({mean_ccd_temp} - robophot_frames.ccd_temp)::numeric) < {delta_temp}")
            '''
            - Оставляешь только те кадры, где CCD-TEMP
            отличается по модулю от mean_ccd_temp
            не более чем на epsilon_temp
            '''
            flats_for_master = pd.read_sql_query(q_get_flats_for_m, eng)
            q_get_m_dark_to_make_flat = (f"SELECT m_frame_path FROM "
                                         f"robophot_master_frames "
                                         f"WHERE m_frame_type = 'm_Dark' AND "
                                         f"{m_flat['m_x_bin']} = m_x_bin AND {m_flat['m_y_bin']} = m_y_bin AND "
                                         f"'{m_flat['m_camera_sn']}' = m_camera_sn "
                                         f"AND {m_flat['m_exp_time']} = m_exp_time AND "
                                         f"abs(({mean_ccd_temp} - ccd_temp_mean)::numeric) < {delta_temp} AND "
                                         f"('{flats_for_master['date_utc'][0]}' - date_make_utc)::INTERVAL"
                                         f" < '2 DAY'::interval")
            m_darks_for_flat = pd.read_sql_query(q_get_m_dark_to_make_flat, eng)
            '''
            - Если после этого в наборе осталось меньше 
            min_images_num (у меня это 5 шт),
            то создание мастер-дарка провалено
            '''
            if flats_for_master.shape[0] < 5 or m_darks_for_flat.shape[0] == 0:
                print("Too few frames to make master flat")
                with eng.begin() as conn:     # TRANSACTION
                    conn.execute(text(set_m_trouble))
                continue
            date_obs = Time(flats_for_master['date_utc'][0])
            date_made = date_obs.to_value(format='fits')
            # path = FLAT_ROOT + f'{m_flat["m_exp_time"]}_{mean_ccd_temp}.fits.gz'
            path_root = rf'D:\{date_obs.to_value("ymdhms")[0]}\{date_obs.to_value("iso", subfmt="date")}\UNKNOWN'
            if not os.path.exists(path_root):
                os.makedirs(path_root)
            path = rf'{path_root}\MFLAT_{date_made.replace(":", "-")}_{m_flat["m_frame_filter"]}_{m_flat["m_exp_time"]}_{mean_ccd_temp}.fits.gz'
            for d_indx, dark in m_darks_for_flat.iterrows():
                if make_master_flat(flats_for_master['frame_path'], dark['m_frame_path'], mean_ccd_temp, 
                                    Time.now().to_value(format='fits'), path):
                    q_update_m_flat = ("UPDATE robophot_master_frames SET "
                                       "(m_frame_path, ccd_temp_mean, date_make_utc) = "
                                       f"('{path}'::text, '{mean_ccd_temp}'::numeric, "
                                       f"'{date_made}'::timestamp) WHERE "
                                       f"m_frame_id = {m_flat['m_frame_id']}")
                    with eng.begin() as conn:     # TRANSACTION
                        conn.execute(text(q_update_m_flat))
                    # con.execute(q_update_m_flat)
                    # eng.execute(text(q_update_m_flat))
                    print(f"Made master flat {path}, used dark {dark['m_frame_path']}")
                    break
                else:
                    print("Failed while make master flat, trying next dark or exit")
                    with eng.begin() as conn:     # TRANSACTION
                        conn.execute(text(set_m_trouble))
        else:
            print('No frames to calibrate')
    eng.dispose()


while True:
    sisyphus()
    time.sleep(10)
