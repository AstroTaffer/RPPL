from rppl_data import RPPLData
from utils import create_empty_configs


create_empty_configs("APM-1.2M.json", "config_calc.json")


"""
H:/STUDY TO/Астрономия/4_НИР/

TODO:

TODO BONUS:
    Добавить файловые логи
    Написасть функцию построения FITS-файла
    Подумать над реализацией функции автосохранения. Возможно, стоит использовать один общий переключатель?

Если я не избавлюсь от учёта количества снимков для построчного считывания промежуточных данных, то
мне придётся вызывать инвентаризацию снимков прежде чем я смогу нормально считать промежуточные данные

# FIXME OLD
#   Plotting functions also have confusing names and bad code
#   ROMS critetia needs to be rewritten and tested on higher quality data
#   SINCE I'M CALCULATING "BJD % 1" I'M STUCK WITHIN BOUNDS OF ONE JULIAN DAY
"""
