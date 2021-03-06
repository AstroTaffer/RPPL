from rppl_data import RPPLData
from utils import plot_results

plot_results()

# da_ob = RPPLData(head_config="APM-1.2M", calc_config="basic_calc")
# da_ob.account_fits_files()
# da_ob.write_images_table("images_table.txt")
# da_ob.make_sdarks()
# da_ob.make_sflats()
# da_ob.write_sdfimg_table("sdfimg_table.txt")
# print(da_ob.sdfimg_table["SIGMA"])

# da_ob.calibrate_images()
# da_ob.write_sdfimg_table("sdfimg_table.txt")


"""
Писать ".json" не нужно, а ".txt" - нужно. Это может путать.


TODO:

TODO BONUS:
    Если в хедере нет фильтра, то в таблицу учета снимков пишется условный фильтр-заглушка
    Добавить файловые логи
    Написасть функцию построения FITS-файла
    Подумать над реализацией функции автосохранения. Возможно, стоит использовать один общий переключатель?
    Если я не избавлюсь от учёта количества снимков для построчного считывания промежуточных данных, ---
    --- то мне придётся вызывать инвентаризацию снимков прежде чем я смогу нормально считать промежуточные данные

# FIXME OLD
#   Plotting functions also have confusing names and bad code
#   ROMS critetia needs to be rewritten and tested on higher quality data
#   Since I'm calculating "BJD % 1" I'm stuck within bounds of one julian day

Бонус
    Визуально проверить снимки на адекватную калибровку
"""
