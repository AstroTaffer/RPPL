from rppl_data import RPPLData

da_ob = RPPLData(head_config="APM-1.2M", calc_config="basic_calc")
da_ob.fits_files_accounting()
da_ob.write_ffa_results()

"""
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
