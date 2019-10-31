import time
import numpy as np
import pandas as pd
import config as cfg
from hexgrid import HexGrid
from helpers import TightBindingHelpers
from lattice_maker import LatticeConstructor
from tight_binding_calculations import TightBinding

class ExecuteTightBindingCallculations(object):

	def __init__(self):
		self.__hexagonal_lattice = HexGrid()
		self.__tight_binding = TightBinding()
		self.__rectangle_lattice = LatticeConstructor()
		self.__configuration = cfg.configuration['parametrization']

	def __call_rectangle_lattice(self, defects: bool = False, number_of_defects: (int, bool) = None, **kwargs) -> (pd.DataFrame, int):
	    lattice = self.__rectangle_lattice.construct_skeleton(**kwargs)
	    if defects == True:
	        lattice = self.__rectangle_lattice.add_defects(lattice, number_of_defects)
	    lattice_data_frame_format = self.__rectangle_lattice.dump_to_dframe(lattice)
	    shape = int(lattice.shape[0])
	    return lattice, lattice_data_frame_format, shape

	def __construct_rectangle_lattice(self):
		distance = self.__configuration['distance']
		vertical_num_of_steps = self.__configuration['vertical_num_of_steps']
		horizontal_num_of_steps = self.__configuration['horizontal_num_of_steps']

		lattice, lattice_df_format, dims = self.__call_rectangle_lattice(a=distance,
				                                                       vertical_num_of_steps=vertical_num_of_steps,
				                                                       horizontal_num_of_steps=horizontal_num_of_steps)
		return lattice, lattice_df_format, dims

	def __call_hexagonal_lattice(self, distance, x_num_of_steps):
	    hexgrid = self.__hexagonal_lattice.generate_hexgrid(distance, x_num_of_steps)
	    df_hex = self.__hexagonal_lattice.dump_to_dframe(hexgrid)
	    shape = int(df_hex.shape[0])
	    return hexgrid, df_hex, shape

	def __construct_hexagonal_lattice(self):
		distance = self.__configuration['distance']
		x_num_of_steps = self.__configuration['x_num_of_steps']
		lattice, lattice_df_format, dims = self.__call_hexagonal_lattice(distance=distance, 
																	   x_num_of_steps=x_num_of_steps)
		return lattice, lattice_df_format, dims

	def choose_type_of_lattice(self):
		lattice_type = self.__configuration['lattice_type']
		if lattice_type == 'hexagonal':
			return self.__construct_hexagonal_lattice()
		else:
			return self.__construct_rectangle_lattice()

	def call_tight_binding_calculation(self, dimension, lattice_df_format):
		tb = TightBinding()
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
		energy, wave_function = tb.calculate_eigenvalues_ang_eigenvectors(dimension=dimension,
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

	def calculate_DOS(self, eigen_energies):
		tb = TightBinding()
		start = self.__configuration['start']
		end = self.__configuration['stop']
		step = self.__configuration['step']
		gauss_sigma = self.__configuration['gauss_sigma']
		E = np.arange(start, end, step)
		density_of_states = tb.evaluate_density_of_states(eigen_energies, E, a=gauss_sigma)
		return density_of_states

def main():
	print('_______________________________TIGHT_BINDING_CALCULATIONS_______________________________\n')
	execute = ExecuteTightBindingCallculations()
	helpers = TightBindingHelpers()

	start = time.time()

	lattice, lattice_df_format, dimension = execute.choose_type_of_lattice()
	energies, wave_function = execute.call_tight_binding_calculation(dimension, lattice_df_format)
	density_of_states = execute.calculate_DOS(energies)

	end = time.time()
	print('Calculation time for ' + str(dimension)  + ' atoms: ', (end - start) / 60., ' minutes')
	file_name = helpers.get_file_name(dimension)
	helpers.create_saving_folder()
	helpers.save_numerical_results(file_name, energies)
	helpers.plot_DOS(file_name, density_of_states)
	helpers.plot_hex_lattice(lattice, file_name)
	return


if __name__ == '__main__':
	exit(main())