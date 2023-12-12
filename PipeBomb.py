import time

import pandas as pd
from astropy.time import Time
import utils
from sqlalchemy import text
from Astrometry import DoAstrometry
from photometric_calibration import apply_dark, apply_flat, make_master_dark, make_master_flat

time_wait_sec = 10
# root = r'D:\RoboPhotData\Images'
# BD_path = root + r'\Calibrated.csv'
DARK_ROOT = r'D:\_Darks\mdark_'   # mflat_exp_temp
FLAT_ROOT = r'D:\_Flats\mflat_'
ID = 'frame_id'
PATH_NAME = 'frame_path'
COL_PATH_NAME = 'calibration_frame_path'
DO_DARK_NAME = 'is_do_dark'
DO_FLAT_NAME = 'is_do_flat'
DO_ASTROMETRY_NAME = 'is_do_astrometry'
TROUBLES_NAME = 'trouble_count'
# names = [PATH_NAME, DO_DARK_NAME, DO_FLAT_NAME, DO_ASTROMETRY_NAME, TROUBLES_NAME]
API_KEY = 'hipfhzhlzygnlvix'
delta_temp = 1
med_fwhm = 0
med_ell = 0
med_bkg = 0
med_zeropoi = 0
sex_path = ''

def sisyphus():
    # timer.cancel()

    print(f'New sisyphus tic at {Time.now()}')
    # read bd

    # if not eng:
    #     print('Error in db connection')
    #     return
    # with (utils.connect_to_db() as eng):
    eng = utils.connect_to_db()
    # con = eng.connect()
    print('Connected to robophot db')
    # query = "SELECT get_frames_for_cal() FROM robophot_frames"
    query = ("SELECT frame_id, frame_path, calibration_frame_path, coord2000, "
             "frame_filter, date_utc, ccd_temp, is_do_dark, is_do_flat, is_do_astrometry, trouble_count "
             "FROM robophot_frames, robophot_tasks "
             "WHERE fk_task_id = task_id AND (frame_type = 'Object' OR frame_type = 'Test') AND "
             "(NOT is_do_dark OR NOT is_do_flat OR NOT is_do_astrometry OR NOT trouble_count > 3) "
             "ORDER BY (CASE "
             "WHEN frame_filter = 'i' THEN 1 "
             "WHEN frame_filter = 'r' THEN 2 "
             "WHEN frame_filter = 'g' THEN 3 END), frame_id")
    table = pd.read_sql_query(query, eng)

    # apply calibration
    if table.shape[0] > 0:
        print(f'Start calibrate, got {table.shape[0]} frames')
        ass = DoAstrometry(API_KEY)
        for indx, row in table.iterrows():
            print("*"*42)
            print(f"Calibration: working on frame #{row[ID]}")
            if not row[DO_DARK_NAME]:
                # q_get_for_do_dark = f"SELECT get_m_dark_frames({row[ID]}, {delta_temp})"
                q_get_for_do_dark = (f"SELECT m_frame_path FROM "
                                     f"robophot_master_frames, robophot_frames, robophot_tasks "
                                     f"WHERE {row[ID]} = frame_id AND robophot_frames.fk_task_id = task_id AND "
                                     f"m_frame_type = 'm_Dark' AND "
                                     f"ABS(ccd_temp-ccd_temp_mean) < {delta_temp} AND camera_sn = m_camera_sn AND "
                                     f"x_bin = m_x_bin AND y_bin = m_y_bin AND exp_time = m_exp_time "
                                     f"ORDER BY ABS(EXTRACT(DAY FROM (date_make_utc - date_utc)::INTERVAL)) "
                                     f"LIMIT 4")
                darks_table = pd.read_sql_query(q_get_for_do_dark, eng)
                for d_index, dark in darks_table.iterrows():
                    if apply_dark(row[PATH_NAME], dark['m_frame_path'], row[COL_PATH_NAME]):
                        row[DO_DARK_NAME] = True
                        print('Dark applied')
                        break
                    else:
                        row[TROUBLES_NAME] += 1
                        print(f"Can\'t apply dark {dark}")

            if row[DO_DARK_NAME] and not row[DO_FLAT_NAME]:
                # q_get_for_do_flat = f"SELECT get_m_flat_frames({row[ID]})"
                q_get_for_do_flat = (
                    f"SELECT m_frame_path FROM robophot_master_frames, robophot_frames, robophot_tasks "
                    f"WHERE {row[ID]} = frame_id AND robophot_frames.fk_task_id = task_id AND "
                    f"m_frame_type = 'm_Flat' AND "
                    f"m_frame_filter = frame_filter AND "
                    f"camera_sn = m_camera_sn AND x_bin = m_x_bin AND y_bin = m_y_bin "
                    f"ORDER BY ABS(EXTRACT(DAY FROM (date_make_utc - date_utc)::INTERVAL))  "
                    f"LIMIT 4")
                flats_table = pd.read_sql_query(q_get_for_do_flat, eng)
                for f_indx, flat in flats_table.iterrows():
                    if apply_flat(row[PATH_NAME], flat['m_frame_path'], row[COL_PATH_NAME]):
                        row[DO_FLAT_NAME] = True
                        print('Flat applied')
                        break
                    else:
                        row[TROUBLES_NAME] += 1
                        print(f"Can\'t apply flat {flat}")

            if row[DO_DARK_NAME] and row[DO_FLAT_NAME] and not row[DO_ASTROMETRY_NAME]:
                if ass.compute(row[COL_PATH_NAME]):
                    row[DO_ASTROMETRY_NAME] = True
                    print('Made flat')
                    med_fwhm, med_ell, med_bkg, med_zeropoi, sex_path = utils.do_sex(row[COL_PATH_NAME])
                else:
                    row[DO_ASTROMETRY_NAME] = False
                    row[TROUBLES_NAME] += 1
                    print("Can\'t make WCS")
            # update
            q_update = (f"UPDATE robophot_frames SET (sex_path, sex_fwhm, sex_ell, sex_background, sex_zeropoint, "
                        f"{DO_DARK_NAME}, {DO_FLAT_NAME}, {DO_ASTROMETRY_NAME}, "
                        f"{TROUBLES_NAME}) = ({sex_path}, {med_fwhm}, {med_ell}, {med_bkg}, {med_zeropoi}, "
                        f"{row[DO_DARK_NAME]}, {row[DO_FLAT_NAME]}, {row[DO_ASTROMETRY_NAME]}, "
                        f"{row[TROUBLES_NAME]}) WHERE {ID} = {row[ID]}")
            with eng.begin() as conn:     # TRANSACTION
                conn.execute(text(q_update))
            # con.execute(q_update)
            # eng.execute(text(q_update))
    else:
        print('No frames to calibrate')

    # make master dark
    q_get_m_dark_frame = ("SELECT m_frame_id, m_frame_type, "
                          "m_camera_sn, m_x_bin, m_y_bin, m_exp_time FROM "
                          "robophot_master_frames WHERE "
                          "date_make_utc IS NULL AND m_frame_type = 'm_Dark'")
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
        q_get_mean_ccd_temp = ("SELECT round(avg(ccd_temp)::numeric, 3) FROM robophot_frames, robophot_tasks WHERE "
                               f" frame_type = 'Dark' AND "
                               f"fk_task_id = task_id AND  "
                               f"{m_dark['m_x_bin']} = x_bin AND {m_dark['m_y_bin']} = y_bin AND "
                               f"'{m_dark['m_camera_sn']}' = camera_sn AND {m_dark['m_exp_time']} = exp_time AND "
                               f"(now() - date_utc) < '1 MONTH'::interval")
        mean_ccd_temp = pd.read_sql_query(q_get_mean_ccd_temp, eng)['round'][0]
        q_get_darks_for_m = (f"SELECT frame_path FROM "
                             f"robophot_frames, robophot_tasks "
                             f"WHERE frame_type = 'Dark' AND "
                             f"fk_task_id = task_id AND "
                             f"{m_dark['m_x_bin']} = x_bin AND {m_dark['m_y_bin']} = y_bin AND "
                             f"'{m_dark['m_camera_sn']}' = camera_sn AND {m_dark['m_exp_time']} = exp_time AND "
                             f"abs(({mean_ccd_temp} - robophot_frames.ccd_temp)::numeric) < {delta_temp} AND "
                             f"(now() - date_utc) < '1 MONTH'::interval")
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
            continue
        date_made = Time.now().to_value(format='fits')
        path = DARK_ROOT + f'{m_dark["m_exp_time"]}_{mean_ccd_temp}.fits.gz'
        if make_master_dark(darks_for_master['frame_path'], mean_ccd_temp, date_made, path):
            q_update_m_dark = ("UPDATE robophot_master_frames SET (m_frame_path, ccd_temp_mean, date_make_utc) = "
                               f"('{path}'::text, {mean_ccd_temp}::numeric, '{date_made}'::timestamp) WHERE "
                               f"m_frame_id = {m_dark['m_frame_id']}")
            with eng.begin() as conn:     # TRANSACTION
                conn.execute(text(q_update_m_dark))
            # con.execute(q_update_m_dark)
            # eng.execute(text(q_update_m_dark))
            print(f"Made master dark {path}")
        else:
            print("Failed while make master dark")

    # make flat
    q_get_m_flat_frame = ("SELECT m_frame_id, m_frame_type, m_frame_filter, "
                          "m_camera_sn, m_x_bin, m_y_bin, m_exp_time FROM robophot_master_frames WHERE "
                          "date_make_utc IS NULL AND m_frame_type = 'm_Flat'")
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
            q_get_mean_ccd_temp = ("SELECT round(avg(ccd_temp)::numeric, 3) FROM robophot_frames, robophot_tasks WHERE "
                                   f"frame_type = 'Flat' AND "
                                   f"'{m_flat['m_frame_filter']}' = frame_filter AND "
                                   f"fk_task_id = task_id AND "
                                   f"{m_flat['m_x_bin']} = x_bin AND {m_flat['m_y_bin']} = y_bin AND "
                                   f"'{m_flat['m_camera_sn']}' = camera_sn AND {m_flat['m_exp_time']} = exp_time AND "
                                   f"(now() - date_utc) < '1 MONTH'::interval")
            mean_ccd_temp = pd.read_sql_query(q_get_mean_ccd_temp, eng)['round'][0]
            q_get_flats_for_m = (f"SELECT frame_path FROM "
                                 f"robophot_frames, robophot_tasks "
                                 f"WHERE frame_type = 'Flat' AND "
                                 f"'{m_flat['m_frame_filter']}' = frame_filter AND "
                                 f"fk_task_id = task_id AND "
                                 f"{m_flat['m_x_bin']} = x_bin AND {m_flat['m_y_bin']} = y_bin AND "
                                 f"'{m_flat['m_camera_sn']}' = camera_sn AND {m_flat['m_exp_time']} = exp_time AND "
                                 f"abs(({mean_ccd_temp} - robophot_frames.ccd_temp)::numeric) < {delta_temp} AND "
                                 f"(now() - date_utc) < '1 MONTH'::interval")
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
                                         f"(now() - date_make_utc) < '1 MONTH'::interval")
            m_darks_for_flat = pd.read_sql_query(q_get_m_dark_to_make_flat, eng)
            '''
            - Если после этого в наборе осталось меньше 
            min_images_num (у меня это 5 шт),
            то создание мастер-дарка провалено
            '''
            if flats_for_master.shape[0] < 5 or m_darks_for_flat.shape[0] == 0:
                print("Too few frames to make master flat")
                continue
            date_made = Time.now().to_value(format='fits')
            path = FLAT_ROOT + f'{m_flat["m_exp_time"]}_{mean_ccd_temp}.fits.gz'
    
            for d_indx, dark in m_darks_for_flat.iterrows():
                if make_master_flat(flats_for_master['frame_path'], dark, mean_ccd_temp, date_made, path):
                    q_update_m_flat = ("UPDATE robophot_master_frames SET "
                                       "(m_frame_path, ccd_temp_mean, date_make_utc) = "
                                       f"('{path}'::text, '{mean_ccd_temp}'::numeric, '{date_made}'::timestamp) WHERE "
                                       f"m_frame_id = {m_flat['m_frame_id']}")
                    with eng.begin() as conn:     # TRANSACTION
                        conn.execute(text(q_update_m_flat))
                    # con.execute(q_update_m_flat)
                    # eng.execute(text(q_update_m_flat))
                    print(f"Made master flat {path}, used dark {dark['frame_path']}")
                    break
                else:
                    print("Failed while make master flat, trying next dark or exit")
        else:
            print('No frames to calibrate')
    eng.dispose()
    # timer.start()


# timer = PerpetualTimer(time_wait_sec, sisyphus)
while True:
    sisyphus()
    time.sleep(60)
