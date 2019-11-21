import os
import datetime
import numpy as np
import config as cfg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class TightBindingHelpers(object):
    """
    Class with helpers method such as plotting and saving results functions of calculations. This method is used in
        ExecuteTightBindingCalculations class.
    """

    def __init__(self, parametrization)->None:
        """
        Method calls configuration file (named config.py)
        """
        self.__configuration = cfg.configuration[parametrization]

    def create_saving_folder(self)->None:
        """
        Method creates (if necessary) folder / directory where results will be saved

        Returns: None

        """
        directory = self.__configuration['saving_directory']
        if not os.path.exists(directory):
            os.makedirs(directory)
        return

    def plot_DOS(self, save_file: str, num_of_atoms: str,density_of_states: np.array)->None:
        """
        Method plots density of states and saves it (as png file).
        Args:
            save_file: name of png file
            density_of_states: array with numerical values of density of states
			num_of_atoms: number of atoms in lattice.

        Returns: None

        """
        start = self.__configuration['start']
        end = self.__configuration['stop']
        step = self.__configuration['step']
        E = np.arange(start, end, step)
        plt.figure(figsize=(13.66, 7.68))
        plt.plot(E, density_of_states)
        plt.axhline(y=0, color='r', linestyle='-')
        plt.xlabel('Energy [a.u.]')
        plt.ylabel('DOS')
        plt.title('Density of states for ' + str(num_of_atoms) + ' atoms')
        plt.savefig(self.__configuration['saving_directory'] + '/__DOS__' + save_file + '.png', dpi=400)
        plt.close()
        return

    def plot_lattice(self, lattice: np.array, save_file:str, num_of_atoms: str, projection: str=None)->None:
        """
        Method plots hexagonal grid.
        Args:
            lattice: array of coordinates of each point of lattice ([[x1,y1,z1], [x2,y2,z2], [...]])
            save_file: name of file
            projection: argument which tells if plot has to be plotted in 3d
            num_of_atoms: number of atoms in lattice.

        Returns: None

        """
        x = []
        y = []
        z = []
        for i in lattice:
            x.append(i[0])
            y.append(i[1])
            z.append([i[2]])
        if projection == '3d':
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(x, y, z, marker='o')
            plt.savefig(self.__configuration['saving_directory'] + '/_lattice_' + save_file + '.png', dpi=200)
            plt.title('Lattice scheme for ' + str(num_of_atoms) + ' atoms')
            plt.axis('off')
            plt.close()
        else:
            plt.figure(figsize=(15, 10))
            plt.scatter(x, y, marker='o')
            plt.savefig(self.__configuration['saving_directory'] + '/_lattice_' + save_file + '.png', dpi=200)
            plt.axis('off')
            plt.close()
        return

    def save_numerical_results(self, save_file: str, eigen_energies: np.array)->None:
        """
        Method saves numerical results - eigen energies - into txt file
        Args:
            save_file: name of file
            eigen_energies: array of eigen energies calculated by diagonalization of interaction matrix.

        Returns: None

        """
        with open(self.__configuration['saving_directory'] + '/' + save_file + '.txt', "w") as file:
            for eigen_energy in eigen_energies:
                file.write(str(eigen_energy) + "\n")
        return

    def get_file_name(self, num_of_atoms)->str:
        """
        Method creates name of saving file (txt and png) basing on chosen parameters.
        Args:
            num_of_atoms: number of atoms in lattice.

        Returns: name of saving file

        """
        save_file = None
        sigma = self.__configuration['sigma']
        gauss_sigma = self.__configuration['gauss_sigma']
        lattice_type = self.__configuration['lattice_type']
        x_num_of_steps = self.__configuration['x_num_of_steps']
        calculation_type = self.__configuration['calculation_type']

        if x_num_of_steps != None:
            x_num_of_steps = self.__configuration['x_num_of_steps']
            save_file = str(datetime.datetime.now()) + '_' + calculation_type + \
                        '_sigma=' + \
                        str(sigma) + \
                        '_gauss_sigma=' + \
                         str(gauss_sigma) + \
                        '_num_of_atoms=' + \
                        str(num_of_atoms) + '_' +\
                        str(lattice_type)
        elif x_num_of_steps == None:
            vertical_num_of_steps = self.__configuration['vertical_num_of_steps']
            horizontal_num_of_steps = self.__configuration['horizontal_num_of_steps']

            save_file = str(datetime.datetime.now()) + '_' + calculation_type + \
                        '_sigma=' + \
                         str(sigma) + \
                         '_gauss_sigma=' + \
                         str(gauss_sigma) + \
                         '_width=' + \
                          str(vertical_num_of_steps) + \
                          '_length=' + \
                          str(horizontal_num_of_steps) + \
                          '_num_of_atoms=' + \
                          str(num_of_atoms) + '_' +\
                        str(lattice_type)
        return save_file

