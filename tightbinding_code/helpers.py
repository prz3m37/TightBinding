import os
import glob
import datetime
import numpy as np
import config as cfg
import matplotlib.pyplot as plt


class TightBindingHelpers(object):
    """
    Class with helpers method such as plotting and saving results functions of calculations. This method is used in
        ExecuteTightBindingCalculations class.
    """

    def __init__(self) -> None:
        """
        Method calls configuration file (named config.py)
        """
        self.__log_file = None
        self.__directory = "./tightbinding_results"

    def create_saving_folder(self) -> None:
        """
        Method creates (if necessary) folder / directory where results will be saved

        Returns: None

        """
        if not os.path.exists(self.__directory):
            os.makedirs(self.__directory)
        return

    def create_logfile(self):
        self.__log_file = open("log_file.txt", "a+")
        return

    def close_logfile(self):
        self.__log_file.close()
        return

    def save_log(self, message):
        msg = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S") + " " + \
                  + message
        self.__log_file.write(msg)
        return

    @staticmethod
    def generate_id_key(title):
        time = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        saving_key = str(time) + " " + str(title)
        return saving_key

    def save_all_params_to_file(self, title, configuration):
        saving_key = self.generate_id_key(title)
        params_file = open(saving_key + "_parametrization.txt")
        for key in configuration:
            parameter = configuration[key]
            if type(parameter) == dict:
                params_file.write(str(parameter) + ": ")
                for p in parameter:
                    params_file.write('\t' + "atom type: " + str(p))
                    for k in parameter[p]:
                        params_file.write('\t ' + '\t' + str(k) + ": " + str(parameter[p][k]))
            else:
                params_file.write(key + ": " + str(parameter))
        params_file.close()
        return

    def search_data_on_disc(self, file_name):
        seeking_file = glob.glob('%s/%s.*' % (self.__directory, file_name), recursive=True)
        return seeking_file

    @staticmethod
    def split_close_friends(close_friends, number_of_cpus):
        number_of_close_friends = len(close_friends)
        avg = number_of_close_friends / float(number_of_cpus)
        divided_close_friends = []
        last = 0.0

        while last < number_of_close_friends:
            divided_close_friends.append(close_friends[int(last):int(last + avg)])
            last += avg

        return divided_close_friends
