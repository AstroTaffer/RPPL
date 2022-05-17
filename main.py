from rppl_data import RPPLData


da_ob = RPPLData(images_directory="robophot pics/g/",
                 images_filter="g")
da_ob.pars["config_file"] = "test_config.json"
da_ob.write_config()


"""
TODO:

TODO BONUS:
    move select_fits_file to Utils
    def plot_fits_file    
"""

# FIXME
#   Plotting functions also have confusing names and bad code
#   Current logs are outright awful
#   ROMS critetia needs to be rewritten and tested on higher quality data
#   SINCE I'M CALCULATING "BJD % 1" I'M STUCK WITHIN BOUNDS OF ONE JULIAN DAY
