import time
import numpy as np
import pandas as pd
import config as cfg
from hexgrid import HexGrid
from helpers import TightBindingHelpers
from lattice_maker import LatticeConstructor
from tight_binding_calculations import TightBinding


class ExecuteTightBindingCalculations(object):

    """
    ExecuteTightBindingCalculations is main class where all other classes are called for calculations desired system

    """

    def __init__(self)->None:
        """
        Method calls all necessary classes : TightBinding, LatticeConstructor, HexGrid and LatticeConstructor class.
            The last two classes are unnecessary if user will define his lattice (in required format) by himself.
        """
        self.__hexagonal_lattice = HexGrid()
        self.__tight_binding = TightBinding()
        self.__rectangle_lattice = LatticeConstructor()
        self.__configuration = cfg.configuration['parametrization']

    def __call_rectangle_lattice(self, defects: bool = False, number_of_defects: (int, bool) = None, **kwargs)\
            -> (pd.DataFrame, int):
        """
        Method calls method responsible for creating rectangular lattice (LatticeConstructor)
        Args:
            defects: bool if user wants to add defects int lattice
            number_of_defects: number of atoms to delete
            **kwargs: rest of parameters required to LatticeConstructor methods

        Returns: lattice (array), lattice_df_format (lattice saved in DataFrame format), dims (dimension of lattice
         - number of atoms in lattice).
        """
        lattice = self.__rectangle_lattice.construct_skeleton(**kwargs)
        if defects == True:
            lattice = self.__rectangle_lattice.add_defects(lattice, number_of_defects)
        lattice_data_frame_format = self.__rectangle_lattice.dump_to_dframe(lattice)
        shape = int(lattice.shape[0])
        return lattice, lattice_data_frame_format, shape

    def __construct_rectangle_lattice(self)->(list, pd.DataFrame, int):
        """
        Method responsible for constructing lattice. All required parameters user has to define in configuration file
        Returns: lattice (array), lattice_df_format (lattice saved in DataFrame format), dims (dimension of lattice
         - number of atoms in lattice).

        """
        defects = self.__configuration['defects']
        distance = self.__configuration['distance']
        number_of_defects = self.__configuration['number_of_defects']
        vertical_num_of_steps = self.__configuration['vertical_num_of_steps']
        horizontal_num_of_steps = self.__configuration['horizontal_num_of_steps']

        lattice, lattice_df_format, dims = self.__call_rectangle_lattice(a=distance,
                                                                         vertical_num_of_steps=vertical_num_of_steps,
                                                                         horizontal_num_of_steps=horizontal_num_of_steps,
                                                                         defects=defects,
                                                                         number_of_defects=number_of_defects)
        return lattice, lattice_df_format, dims

    def __call_hexagonal_lattice(self, distance: float, x_num_of_steps: int)->(list, pd.DataFrame, int):
        """
        Method calls method responsible for creating hexagonal lattice (HexGrid)
        Args:
            distance: defined distance between atoms
            x_num_of_steps: number of steps to create lattice

        Returns:  hexgrid - lattice (array), df_hex (lattice saved in DataFrame format),
         shape (dimension of lattice - number of atoms in lattice).

        """

        hexgrid = self.__hexagonal_lattice.generate_hexgrid(distance, x_num_of_steps)
        df_hex = self.__hexagonal_lattice.dump_to_dframe(hexgrid)
        shape = int(df_hex.shape[0])
        return hexgrid, df_hex, shape

    def __construct_hexagonal_lattice(self)->(list, pd.DataFrame, int):
        """
        Method responsible for constructing lattice. All required parameters user has to define in configuration file

        Returns: lattice (array), lattice_df_format (lattice saved in DataFrame format), dims (dimension of lattice
         - number of atoms in lattice).

        """
        distance = self.__configuration['distance']
        x_num_of_steps = self.__configuration['x_num_of_steps']
        lattice, lattice_df_format, dims = self.__call_hexagonal_lattice(distance=distance,
                                                                         x_num_of_steps=x_num_of_steps)
        return lattice, lattice_df_format, dims

    def choose_type_of_lattice(self):
        """
        Method responsible for running appropriate lattice constructor class.
        Returns: required lattice form

        """
        lattice_type = self.__configuration['lattice_type']
        if lattice_type == 'hexagonal':
            return self.__construct_hexagonal_lattice()
        else:
            return self.__construct_rectangle_lattice()

    def call_tight_binding_calculation(self, dimension:int, lattice_df_format:pd.DataFrame)->(np.array, np.array):
        """
        Method call Tight Binding class responsible for  calculation of eigen energies and eigen states (vectors) of
        defined lattice.
        Args:
            dimension: number of atoms
            lattice_df_format: lattice in DataFrame format

        Returns: eigen energies and eigen states (vectors)

        """
        sigma = self.__configuration['sigma']
        which = self.__configuration['magnitude']
        distance = self.__configuration['distance']
        atom_store = self.__configuration['diagonal_energies']
        lanczos_vectors = self.__configuration['lanczos_vectors']
        calculation_type = self.__configuration['calculation_type']
        method = self.__configuration['neighbour_calculation_method']
        number_of_friends = self.__configuration['number_of_friends']
        num_of_eigenvalues = self.__configuration['number_of_eigenvalues']
        constants_of_pairs = self.__configuration['interactive_constans']
        energy, wave_function = self.__tight_binding.calculate_eigenvalues_ang_eigenvectors(
                                                                          dimension=dimension,
                                                                          data=lattice_df_format,
                                                                          calculation_type=calculation_type,
                                                                          method=method,
                                                                          distance=distance,
                                                                          number_of_friends=number_of_friends,
                                                                          num_of_eigenvalues=num_of_eigenvalues,
                                                                          sigma=sigma,
                                                                          constants_of_pairs=constants_of_pairs,
                                                                          atom_store=atom_store,
                                                                          which=which,
                                                                          lanczos_vectors=lanczos_vectors)
        return energy, wave_function

    def calculate_DOS(self, eigen_energies:np.array)->np.array:
        """
        Method calculating density of states
        Args:
            eigen_energies: eigen energies of defined system

        Returns: density of states (array)

        """
        start = self.__configuration['start']
        end = self.__configuration['stop']
        step = self.__configuration['step']
        gauss_sigma = self.__configuration['gauss_sigma']
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
        start = self.__configuration['start']
        end = self.__configuration['stop']
        step = self.__configuration['step']
        gauss_sigma = self.__configuration['gauss_sigma']
        E = np.arange(start, end, step)
        projected_density_of_states = self.__tight_binding.evaluate_projected_density_of_states(eigen_energies,
                                                                                                E,
                                                                                                eigen_vectors,
                                                                                                gauss_sigma)
        return projected_density_of_states


def main()->None:
    """
    Function where user calls all required methods and classes
    Returns: None

    """
    print('_________________________________TIGHT_BINDING_CALCULATIONS_________________________________\n')
    execute = ExecuteTightBindingCalculations()
    helpers = TightBindingHelpers()

    start = time.time()

    lattice, lattice_df_format, dimension = execute.choose_type_of_lattice()
    energies, wave_function = execute.call_tight_binding_calculation(dimension, lattice_df_format)
    density_of_states = execute.calculate_DOS(energies)

    end = time.time()
    print('_______________Calculation time for ' + str(dimension) + ' atoms: ', round((end - start)%60., 3), ' minutes_______________\n')
    file_name = helpers.get_file_name(dimension)
    helpers.create_saving_folder()
    helpers.save_numerical_results(file_name, energies)
    helpers.plot_DOS(file_name, density_of_states)
    helpers.plot_hex_lattice(lattice, file_name)
    return


if __name__ == '__main__':
    exit(main())
