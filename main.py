from rppl_data import RPPLData


my_data = RPPLData()
my_data.list_fits_files()
# my_data.make_master_dark()
# my_data.make_master_flat()
my_data.write_images_table("ImagesTable.txt")
my_data.old_apply_calibration()


"""
TODO:
   Функция добавления кадра в таблицу
"""
