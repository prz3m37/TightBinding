from tight_binding import TightBinding
from data_manager import DataManager
import sys
import warnings
import numpy as np
import pandas as pd
import config as cfg
warnings.filterwarnings("ignore")


class ExecuteTightBindingCalculations(object):

    """
    ExecuteTightBindingCalculations is main class where all other classes are called for calculations desired system

    """

    def __init__(self, helpers, configuration, settings) -> None:
        """
        Method calls all necessary classes : TightBinding, LatticeConstructor, HexGrid and LatticeConstructor class.
            The last two classes are unnecessary if user will define his lattice (in required format) by himself.
        """

        self.__helpers = helpers
        self.__settings = settings
        self.__data_manager = DataManager(self.__helpers, self.__settings)
        self.__tight_binding = TightBinding(self.__helpers, configuration)

    def __call_tight_binding_calculation(self, lattice_df_format: pd.DataFrame) -> (np.array, np.array):
        """
        Method call Tight Binding class responsible for  calculation of eigen energies and eigen states (vectors) of
        defined lattice.
        Args:
            lattice_df_format: lattice in DataFrame format

        Returns: eigen energies and eigen states (vectors)
        """
        lattice = lattice_df_format
        number_of_atoms = lattice.shape[0]

        energy, wave_function, interaction_matrix = self.__tight_binding.calculate_eigenvalues_and_eigenvectors(
            data=lattice, number_of_atoms=number_of_atoms,)
        return energy, wave_function, interaction_matrix

    def __calculate_DOS(self, eigen_energies: np.array) -> np.array:
        """
        Method calculating density of states
        Args:
            eigen_energies: eigen energies of defined system

        Returns: density of states (array)

        """
        start = self.__settings['start']
        end = self.__settings['stop']
        step = self.__settings['step']
        gauss_sigma = self.__settings['gauss_sigma']
        E = np.arange(start, end, step)
        density_of_states = self.__tight_binding.evaluate_density_of_states(
            eigen_energies, E, gauss_sigma)
        return density_of_states

    def __calculate_projected_DOS(self, eigen_energies: np.array, eigen_vectors: np.array) -> np.array:
        """
        Method calculating projected density of states
        Args:
            eigen_energies: eigen energies of defined system
            eigen_vectors: eigen states of defined system

        Returns: projected density of states (array)

        """
        start = self.__settings['start']
        end = self.__settings['stop']
        step = self.__settings['step']
        gauss_sigma = self.__settings['gauss_sigma']
        E = np.arange(start, end, step)
        projected_density_of_states = self.__tight_binding.evaluate_projected_density_of_states(eigen_energies,
                                                                                                E,
                                                                                                eigen_vectors,
                                                                                                gauss_sigma)
        return projected_density_of_states

    def execute_tb_calculations(self, lattice):
        try:
            energies, wave_functions, interaction_matrix = self.__call_tight_binding_calculation(
                lattice)
            density_of_states = self.__calculate_DOS(energies)
            projected_density_of_states = self.__calculate_projected_DOS(
                energies, wave_functions)
            return energies, wave_functions, interaction_matrix, density_of_states, projected_density_of_states
        except RuntimeError:
            self.__helpers.save_log("[ERROR]: Factor is exactly singular \n")
            self.__helpers.save_log("[INFO]: Calculations have stopped \n")
            self.__helpers.close_logfile()
            return
        except TypeError:
            self.__helpers.save_log(
                "[ERROR]: No data to calculate. Please check your configuration or input \n")
            self.__helpers.save_log("[INFO]: Calculations have stopped \n")
            self.__helpers.close_logfile()
            return
