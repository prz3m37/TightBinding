import sys
import numpy as np
import pandas as pd
import config as cfg
from tight_binding import TightBinding
from helpers import TightBindingHelpers


class ExecuteTightBindingCalculations():

    """
    ExecuteTightBindingCalculations is main class where all other classes are called for calculations desired system

    """

    def __init__(self)->None:
        """
        Method calls all necessary classes : TightBinding, LatticeConstructor, HexGrid and LatticeConstructor class.
            The last two classes are unnecessary if user will define his lattice (in required format) by himself.
        """
        sys.path.insert(1, '..')
        self.__settings = cfg.settings
        self.parametrization = sys.argv[1]
        self.__tight_binding = TightBinding()
        self.__configuration = cfg.configuration[self.parametrization]


    def call_tight_binding_calculation(self, dimension: int, lattice_df_format: pd.DataFrame)->(np.array, np.array):
        """
        Method call Tight Binding class responsible for  calculation of eigen energies and eigen states (vectors) of
        defined lattice.
        Args:
            dimension: number of atoms
            lattice_df_format: lattice in DataFrame format

        Returns: eigen energies and eigen states (vectors)
        """


        distance = self.__configuration['distance']
        magnitude = self.__configuration['magnitude']
        fermi_level = self.__configuration['fermi_level']
        atom_store = self.__configuration['diagonal_energies']
        lanczos_vectors = self.__configuration['lanczos_vectors']
        calculation_type = self.__configuration['calculation_type']
        number_of_friends = self.__configuration['number_of_friends']
        constants_of_pairs = self.__configuration['interactive_constans']
        num_of_eigenvalues = self.__configuration['number_of_eigenvalues']
        neighbour_calculation_method = self.__configuration['neighbour_calculation_method']
        energy, wave_function = self.__tight_binding.calculate_eigenvalues_ang_eigenvectors(
                                                                          dimension=dimension,
                                                                          data=lattice_df_format,
                                                                          calculation_type=calculation_type,
                                                                          method=neighbour_calculation_method,
                                                                          distance=distance,
                                                                          number_of_friends=number_of_friends,
                                                                          num_of_eigenvalues=num_of_eigenvalues,
                                                                          fermi_level=fermi_level,
                                                                          constants_of_pairs=constants_of_pairs,
                                                                          atom_store=atom_store,
                                                                          magnitude=magnitude,
                                                                          lanczos_vectors=lanczos_vectors)
        return energy, wave_function

    def calculate_DOS(self, eigen_energies:np.array)->np.array:
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
        density_of_states = self.__tight_binding.evaluate_density_of_states(eigen_energies, E, gauss_sigma)
        return density_of_states

    def calculate_projected_DOS(self, eigen_energies:np.array, eigen_vectors:np.array)->np.array:
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

