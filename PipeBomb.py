import threading
import glob
import pandas as pd
import os.path
from Astrometry import DoAstrometry


time_wait_sec = 10
root = r'D:\RoboPhotData\Images'
BD_path = root + r'\Calibrated.csv'
PATH_NAME = 'Path'
DO_DARK_NAME = 'IsDoDark'
DO_FLAT_NAME = 'IsDoFlat'
DO_ASTROMETRY_NAME = 'IsDoAstrometry'
TROUBLES_NAME = 'TroubleCount'
names = [PATH_NAME, DO_DARK_NAME, DO_FLAT_NAME, DO_ASTROMETRY_NAME, TROUBLES_NAME]
API_KEY = 'hipfhzhlzygnlvix'


def Sisyphus():
    timer.cancel()
    # read bd
    if os.path.exists(BD_path):
        table = pd.read_csv(BD_path)
    else:
        table = pd.DataFrame(columns=[names])
    paths_i = glob.glob(root+r'\*\RAW\i\*.fits')
    paths_r = glob.glob(root+r'\*\RAW\r\*.fits')
    paths_g = glob.glob(root+r'\*\RAW\g\*.fits')
    filters = [paths_i, paths_r, paths_g]

    ass = DoAstrometry(API_KEY)
    for fil in filters:
        for file in fil:
            error = 0
            output_path = file.replace('RAW', 'CLB')
            mask = table[PATH_NAME].isin(file)
            row = table[mask]
            if not row.empty():
                if row[TROUBLES_NAME] > 0:
                    print("*"*42)
                    print(f"Calibration: working on {file}")
                    if not row[DO_DARK_NAME]:
                        # DoDark
                        pass
                    if not row[DO_FLAT_NAME]:
                        # DoFlat
                        pass
                    if not row[DO_ASTROMETRY_NAME]:
                        if ass.compute(output_path):
                            row[DO_ASTROMETRY_NAME] = True
                        else:
                            row[DO_ASTROMETRY_NAME] = False
                            error += 1
                    row[TROUBLES_NAME] = error
                    # update
            else:
                print("*"*42)
                print(f"Calibration: working on {file}")
                row = pd.Series(name=names)
                # DoDark
                # DoFlat
                if ass.compute(output_path):
                    row[DO_ASTROMETRY_NAME] = True
                else:
                    row[DO_ASTROMETRY_NAME] = False
                    error += 1
                row[TROUBLES_NAME] = error
                table.add(row)
    table.to_csv(BD_path)
    timer.start()


timer = threading.Timer(time_wait_sec, Sisyphus)
timer.start()
