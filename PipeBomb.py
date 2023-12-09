import threading

import pandas as pd
from astropy.time import Time
import utils
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


def sisyphus():
    timer.cancel()
    # read bd

    # if not eng:
    #     print('Error in db connection')
    #     return
    with (utils.connect_to_db() as eng):
        query = "SELECT get_frames_for_cal() FROM robophot_frames"
        table = pd.read_sql_query(query, eng)

        # apply calibration
        if table.shape[0] > 0:
            ass = DoAstrometry(API_KEY)
            for row in table:
                print("*"*42)
                print(f"Calibration: working on frame #{row[ID]}")
                if not row[DO_DARK_NAME]:
                    q_get_for_do_dark = f"SELECT get_m_dark_frames({row[ID]}, {delta_temp})"
                    darks_table = pd.read_sql_query(q_get_for_do_dark, eng)
                    for dark in darks_table:
                        if apply_dark(row[PATH_NAME], dark, row[COL_PATH_NAME]):
                            row[DO_DARK_NAME] = True
                            print('Dark applied')
                            break
                        else:
                            row[TROUBLES_NAME] += 1
                            print(f"Can\'t apply dark {dark}")

                if row[DO_DARK_NAME] and not row[DO_FLAT_NAME]:
                    q_get_for_do_flat = f"SELECT get_m_flat_frames({row[ID]})"
                    flats_table = pd.read_sql_query(q_get_for_do_flat, eng)
                    for flat in flats_table:
                        if apply_flat(row[PATH_NAME], flat, row[COL_PATH_NAME]):
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
                        # do_sex
                    else:
                        row[DO_ASTROMETRY_NAME] = False
                        row[TROUBLES_NAME] += 1
                        print("Can\'t make WCS")
                # update
                q_update = (f"UPDATE robophot_frames SET {DO_DARK_NAME}, {DO_FLAT_NAME}, {DO_ASTROMETRY_NAME}, "
                            f"{TROUBLES_NAME} = {row[DO_DARK_NAME]}, {row[DO_FLAT_NAME]}, {row[DO_ASTROMETRY_NAME]}, "
                            f"{row[TROUBLES_NAME]} WHERE {ID} = {row[ID]}")

                eng.execute(q_update)

        # make master dark
        q_get_m_dark_frame = ("SELECT m_frame_id, m_frame_type, "
                              "m_camera_sn, m_x_bin, m_y_bin, m_exp_time FROM "
                              "robophot_master_frames WHERE "
                              "date_make_utc IS NULL AND m_frame_type = 'm_Dark'")
        m_darks = pd.read_sql_query(q_get_m_dark_frame, eng)

        '''
        - У кадров должны быть одинаковыми:
          -- IMAGETYP
          -- XBINNING
          -- YBINNING
          -- SERNUM
          -- EXPTIME
        '''
        for m_dark in m_darks:
            q_get_mean_ccd_temp = ("SELECT round(avg(ccd_temp)::numeric, 3) FROM robophot_frames WHERE "
                                   " frame_type = {m_dark['m_frame_type']} AND "
                                   f"fk_task_id = task_id AND "
                                   f"{m_dark['m_x_bin']} = x_bin AND {m_dark['m_y_bin']} = y_bin AND "
                                   f"{m_dark['m_camera_sn']} = camera_sn AND {m_dark['m_exp_time']} = exp_time AND "
                                   f"(now() - date_utc) < '1 MONTH'::interval")
            mean_ccd_temp = pd.read_sql_query(q_get_mean_ccd_temp, eng)['round'][0]
            q_get_darks_for_m = (f"SELECT frame_path FROM "
                                 f"robophot_frames, robophot_tasks "
                                 f"WHERE frame_type = {m_dark['m_frame_type']} AND "
                                 f"fk_task_id = task_id AND "
                                 f"{m_dark['m_x_bin']} = x_bin AND {m_dark['m_y_bin']} = y_bin AND "
                                 f"{m_dark['m_camera_sn']} = camera_sn AND {m_dark['m_exp_time']} = exp_time AND "
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
            if darks_for_master.count() < 5:
                continue
            date_made = Time.now().to_value(format='fits')
            path = DARK_ROOT + f'{m_dark["m_exp_time"]}_{mean_ccd_temp}.fits.gz'
            if make_master_dark(darks_for_master['frame_path'], mean_ccd_temp, date_made, path):
                q_update_m_dark = ("UPDATE robophot_master_frames SET (m_frame_path, ccd_temp_mean, date_make_utc) = "
                                   f"({path}::text, {mean_ccd_temp}::numeric, {date_made}::timestamp) WHERE "
                                   f"m_frame_id = {m_dark['m_frame_id']}")
                eng.execute(q_update_m_dark)
            else:
                pass

        # make flat
        q_get_m_flat_frame = ("SELECT m_frame_id, m_frame_type, m_frame_filter, "
                              "m_camera_sn, m_x_bin, m_y_bin, m_exp_time FROM robophot_master_frames WHERE "
                              "date_make_utc IS NULL AND m_frame_type = 'm_Flat'")
        m_flats = pd.read_sql_query(q_get_m_flat_frame, eng)

        '''
        - У кадров должны быть одинаковыми:
          -- IMAGETYP
          -- XBINNING
          -- YBINNING
          -- SERNUM
          -- EXPTIME
        '''
        for m_flat in m_flats:
            q_get_mean_ccd_temp = ("SELECT round(avg(ccd_temp)::numeric, 3) FROM robophot_frames WHERE "
                                   " frame_type = {m_flat['m_frame_type']} AND "
                                   "{m_flat['m_frame_filter']} = frame_filter AND "
                                   f"fk_task_id = task_id AND "
                                   f"{m_flat['m_x_bin']} = x_bin AND {m_flat['m_y_bin']} = y_bin AND "
                                   f"{m_flat['m_camera_sn']} = camera_sn AND {m_flat['m_exp_time']} = exp_time AND "
                                   f"(now() - date_utc) < '1 MONTH'::interval")
            mean_ccd_temp = pd.read_sql_query(q_get_mean_ccd_temp, eng)['round'][0]
            q_get_flats_for_m = (f"SELECT frame_path FROM "
                                 f"robophot_frames, robophot_tasks "
                                 f"WHERE frame_type = {m_flat['m_frame_type']} AND "
                                 f"{m_flat['m_frame_filter']} = frame_filter AND "
                                 f"fk_task_id = task_id AND "
                                 f"{m_flat['m_x_bin']} = x_bin AND {m_flat['m_y_bin']} = y_bin AND "
                                 f"{m_flat['m_camera_sn']} = camera_sn AND {m_flat['m_exp_time']} = exp_time AND "
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
                                         f"{m_flat['m_camera_sn']} = m_camera_sn "
                                         f"AND {m_flat['m_exp_time']} = m_exp_time AND "
                                         f"abs(({mean_ccd_temp} - ccd_temp_mean)::numeric) < {delta_temp} AND "
                                         f"(now() - date_make_utc) < '1 MONTH'::interval")
            m_darks_for_flat = pd.read_sql_query(q_get_m_dark_to_make_flat, eng)
            '''
            - Если после этого в наборе осталось меньше 
            min_images_num (у меня это 5 шт),
            то создание мастер-дарка провалено
            '''
            if flats_for_master.count() < 5 or m_darks_for_flat.count() == 0:
                continue
            date_made = Time.now().to_value(format='fits')
            path = FLAT_ROOT + f'{m_flat["m_exp_time"]}_{mean_ccd_temp}.fits.gz'

            for dark in m_darks_for_flat:
                if make_master_flat(flats_for_master['frame_path'], dark, mean_ccd_temp, date_made, path):
                    q_update_m_dark = ("UPDATE robophot_master_frames SET "
                                       "(m_frame_path, ccd_temp_mean, date_make_utc) = "
                                       f"({path}::text, {mean_ccd_temp}::numeric, {date_made}::timestamp) WHERE "
                                       f"m_frame_id = {m_flat['m_frame_id']}")
                    eng.execute(q_update_m_dark)
                    break
                else:
                    pass
    timer.start()


timer = threading.Timer(time_wait_sec, sisyphus)
timer.start()
