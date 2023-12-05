import threading
import glob
import pandas as pd
import os.path

import utils
from Astrometry import DoAstrometry


time_wait_sec = 10
root = r'D:\RoboPhotData\Images'
BD_path = root + r'\Calibrated.csv'
ID = 'frame_id'
PATH_NAME = 'frame_path'
COL_PATH_NAME = 'calibration_frame_path'
DO_DARK_NAME = 'is_do_dark'
DO_FLAT_NAME = 'is_do_flat'
DO_ASTROMETRY_NAME = 'is_do_astrometry'
TROUBLES_NAME = 'trouble_count'
names = [PATH_NAME, DO_DARK_NAME, DO_FLAT_NAME, DO_ASTROMETRY_NAME, TROUBLES_NAME]
API_KEY = 'hipfhzhlzygnlvix'


def Sisyphus():
    timer.cancel()
    # read bd
    if not con:
        print('Error in db connection')
        return
    query = "SELECT get_frames_for_cal() FROM robophot_frames"
    table = pd.read_sql_query(query, con)
    if table.shape[0] == 0:
        return

    ass = DoAstrometry(API_KEY)
    for row in table:
        print("*"*42)
        print(f"Calibration: working on frame #{row[ID]}")
        if row[TROUBLES_NAME] > 3:
            print(f"Too many troubles, next frame")
            continue
        error = 0
        if not row[DO_DARK_NAME]:
            # DoDark
            if not row[DO_FLAT_NAME]:
                # DoFlat
                if not row[DO_ASTROMETRY_NAME]:
                    if ass.compute(row[COL_PATH_NAME]):
                        row[DO_ASTROMETRY_NAME] = True
                        # do_sex
                    else:
                        row[DO_ASTROMETRY_NAME] = False
                        error += 1
        row[TROUBLES_NAME] = error
        # update

    timer.start()


timer = threading.Timer(time_wait_sec, Sisyphus)
con = utils.connect_to_db()
con.autocommit = True
timer.start()
