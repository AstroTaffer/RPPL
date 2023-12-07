import threading
import pandas as pd

import utils
from Astrometry import DoAstrometry
from photometric_calibration import apply_dark, apply_flat

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
    with utils.connect_to_db() as eng:
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
                q_update = (f"UPDATE robophpt_frames SET {DO_DARK_NAME}, {DO_FLAT_NAME}, {DO_ASTROMETRY_NAME}, "
                            f"{TROUBLES_NAME} = {row[DO_DARK_NAME]}, {row[DO_FLAT_NAME]}, {row[DO_ASTROMETRY_NAME]}, "
                            f"{row[TROUBLES_NAME]} WHERE ID = {row[ID]}")

                eng.execute(q_update)
        # make master files
        q_get_m_dark_frame = "SELECT "
        m_darks = pd.read_sql_query(q_get_m_dark_frame, eng)
        for m_dark in m_darks:
            q_get_darks_for_m = "SELECT "
            darks_for_master = pd.read_sql_query(q_get_darks_for_m, eng)
            # make_master_dark()
    timer.start()


timer = threading.Timer(time_wait_sec, sisyphus)
timer.start()
