import json


class _RPPLSubIOS:
    def __init__(self):
        self.pars = {}

    def read_config(self):
        with open(self.pars["config_file"], "r") as confile:
            buff_dict = json.load(confile)
            for _ in buff_dict.keys():
                if _ in self.pars:
                    self.pars[_] = buff_dict[_]
                else:
                    print(f"Parameter {_} not recognised --- skipped")

    def write_config(self):
        if self.pars["config_file"] is None:
            print('Parameter config_file not set --- set "config.json" by default')
            self.pars["config_file"] = "config.json"
        with open(self.pars["config_file"], "w") as confile:
            json.dump(self.pars, confile, indent=4)


"""
read_app_results
write_app_results

result file name = f"[RESULT_NAME]_{config_file[:-5]}.txt"
"""