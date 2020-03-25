import pickle
import numpy as np
from typing import TextIO
from tight_binding_base_connection import DataBaseConnection


class DataManager(object):

    def __init__(self, helpers, settings):
        self.__helpers = helpers
        self.__settings = settings
        self.__dbconnection = DataBaseConnection(settings, helpers)

    def __save_results_into_database(self, data):
        table = self.__settings["db table"]
        data_identification = self.__settings["data_identification"]
        save_query = """INSERT INTO {} (Lattice identification,
                                        Eigen energy, 
                                        Wave function, 
                                        Slater-Koster matrix, 
                                        DOS, 
                                        Projected DOS, 
                                        Configuration) 
                                VALUES ({}, %s, %s, %s, %s, %s, %s) """.format(table, data_identification)
        query = (save_query, data)
        return query

    def __download_required_data(self):
        select = self.__settings["select"]
        table = self.__settings["db table"]
        identifier = self.__helpers["identifier"]
        if identifier is None:
            load_query = """SELECT {} from {}""".format(select, table)
        else:
            load_query = """SELECT {} from {} WHERE {}""".format(
                select, table, identifier)
        return load_query

    @staticmethod
    def __save_as_pickle(matrix):
        with open("interaction_matrix", 'wb') as outfile:
            pickled_matrix = pickle.dump(
                matrix, outfile, pickle.HIGHEST_PROTOCOL)
        return pickled_matrix

    @staticmethod
    def __load_as_pickle(pickled_matrix):
        with open(pickled_matrix, 'rb') as infile:
            matrix = pickle.load(infile)
        return matrix

    def __save_numerical_results(self, title, eigen_energies: np.array) -> TextIO:
        """
        Method saves numerical results - eigen energies - into txt file
        Args:
            title: name of file
            eigen_energies: array of eigen energies calculated by diagonalization of interaction matrix.

        Returns: None

        """
        saving_key = self.__helpers.generate_id_key(title)
        with open(self.__helpers.__directory + "/" + saving_key + '_eigenvalues.txt', "w") as file:
            for eigen_energy in eigen_energies:
                file.write(str(eigen_energy) + "\n")
        return file

    def __save_data_locally(self, energies, wave_functions, interaction_matrix, density_of_states,
                            projected_density_of_states, configuration):

        energy_file = self.__save_numerical_results("eigen_energies", energies)
        wave_functions_file = self.__save_as_pickle(wave_functions)
        interaction_matrix_file = self.__save_as_pickle(interaction_matrix)
        dos_file = self.__save_numerical_results("DOS", density_of_states)
        p_dos_file = self.__save_numerical_results(
            "PDOS", projected_density_of_states)
        configuration_file = self.__helpers.save_all_params_to_file(
            "parametrization_file", configuration)
        data_to_save = (energy_file, wave_functions_file,
                        interaction_matrix_file, dos_file, p_dos_file, configuration_file)
        return data_to_save

    def save_data(self, energies, wave_functions, interaction_matrix,
                  density_of_states, projected_density_of_states, configuration):
        data_source = self.__settings["data_source"]
        self.__helpers.save_log('[INFO]: Saving results locally \n')
        data_to_save = self.__save_data_locally(energies,
                                                wave_functions,
                                                interaction_matrix,
                                                density_of_states,
                                                projected_density_of_states,
                                                configuration)
        if data_source == "data base":
            self.__helpers.save_log('[INFO]: Saving results into data base\n')
            save_data_query = self.__save_results_into_database(data_to_save)
            self.__dbconnection.execute_query(save_data_query, "save")
        else:
            pass
        return

    def load_data(self):
        self.__helpers.save_log('[INFO]: Loading results from data base\n')
        data_source = self.__settings["data_source"]
        if data_source == "data base":
            data = self.__download_required_data()
            required_data = self.__dbconnection.execute_query(data, "load")
        else:
            seeking_file = self.__settings["data_identification"]
            required_data = self.__helpers.search_data_on_disc(seeking_file)
        return required_data
