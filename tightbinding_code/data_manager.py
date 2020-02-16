import pickle
import numpy as np
from typing import TextIO
from tight_binding_base_connection import DataBaseConnection


class DataManager(object):

    def __init__(self, helpers, settings):
        self.__helpers = helpers
        self.__settings = settings
        self.__dbconnection = DataBaseConnection()

    def __save_results_into_database(self, data):
        query = ""
        return query

    def __download_required_data(self):
        return ""

    @staticmethod
    def __save_as_pickle(matrix):
        with open("interaction_matrix", 'wb') as outfile:
            pickled_matrix = pickle.dump(matrix, outfile, pickle.HIGHEST_PROTOCOL)
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
        saving_key = self.__helpers.generate_saving_key(title)
        with open(self.__helpers.__directory + "/" + saving_key + '_eigenvalues.txt', "w") as file:
            for eigen_energy in eigen_energies:
                file.write(str(eigen_energy) + "\n")
        return file

    def save_data(self, **kwargs):
        data_source = self.__settings["data_source"]
        self.__helpers.save_log('[INFO]: Saving data locally \n')
        interaction_marix = self.__save_as_pickle(**kwargs)
        eigenenergies = self.__save_numerical_results(**kwargs)
        if data_source == "data base":
            self.__helpers.save_log('[INFO]: Saving data into database \n')
            query_for_eigenvalues = self.__save_results_into_database(eigenenergies)
            query_for_interaction_marix = self.__save_results_into_database(interaction_marix)
            self.__dbconnection.execute_query(query_for_eigenvalues)
            self.__dbconnection.execute_query(query_for_interaction_marix)
        return

    def load_data(self, **kwargs):
        data_source = self.__settings["data_source"]
        if data_source == "data base":
            self.__helpers.save_log('[INFO]: Loading data from data base \n')
            query = self.__download_required_data()
            pickled_matrix = self.__dbconnection.execute_query(query)
            interaction_matrix = self.__load_as_pickle(pickled_matrix)
        else:
            self.__helpers.save_log('[INFO]: Loading data from disk \n')
            interaction_matrix = self.__load_as_pickle(**kwargs)
        return interaction_matrix

# TODO: Maybe check if loaded matrix is proper to "use"