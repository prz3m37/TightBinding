import os
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

    # def plot_DOS(self, save_file: str, num_of_atoms: str, density_of_states: np.array) -> None:
    #     """
    #     Method plots density of states and saves it (as png file).
    #     Args:
    #         save_file: name of png file
    #         num_of_atoms: number of atoms in lattice.
    #         density_of_states: array with numerical values of density of states
    #
    #     Returns: None
    #
    #     """
    #     start = self.__settings['start']
    #     end = self.__settings['stop']
    #     step = self.__settings['step']
    #     E = np.arange(start, end, step)
    #     plt.figure(figsize=(13.66, 7.68))
    #     plt.plot(E, density_of_states)
    #     plt.axhline(y=0, color='r', linestyle='-')
    #     plt.xlabel('Energy [a.u.]')
    #     plt.ylabel('DOS')
    #     plt.title('Density of states for ' + str(num_of_atoms) + ' atoms')
    #     plt.savefig(self.__directory + '/__DOS__' + save_file + '.png', dpi=400)
    #     plt.close()
    #     return

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
    def generate_saving_key(title, save_type="directory"):
        time = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        if save_type == "database":
            saving_key = [time, title]
        else:
            saving_key = str(time) + " " + str(title)
        return saving_key

    def save_all_params_to_file(self, parametrization, title):
        saving_key = self.generate_saving_key(title)
        params_file = open(saving_key + "_parametrization.txt")
        for key in parametrization:
            parameter = parametrization[key]
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

