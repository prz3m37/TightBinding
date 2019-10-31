import os
import numpy as np
import config as cfg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class TightBindingHelpers(object):

	def __init__(self):
		self.__configuration = cfg.configuration['parametrization']

	def create_saving_folder(self):
		directory = self.__configuration['saving_directory'] 
		if not os.path.exists(directory):
			os.makedirs(directory)
		return 

	def plot_DOS(self, save_file, density_of_states):
		start = self.__configuration['start']
		end = self.__configuration['stop']
		step = self.__configuration['step']
		gauss_sigma = self.__configuration['gauss_sigma']
		E = np.arange(start, end, step)
		plt.figure(figsize=(25,15))
		plt.plot(E, density_of_states)
		plt.axhline(y=0, color='r', linestyle='-')
		plt.xlabel('Energy')
		plt.ylabel('DOS')
		plt.savefig(self.__configuration['saving_directory']  + '/__DOS__' + save_file + '.png', dpi=200)
		plt.close()
		return

	def plot_hex_lattice(self, lattice, save_file, projection=None):
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
			plt.savefig(self.__configuration['saving_directory']  + '/_lattice_' + save_file + '.png', dpi=200)
			plt.close()
		else:
			plt.figure(figsize=(15,10))
			plt.scatter(x, y, marker='o')
			plt.savefig(self.__configuration['saving_directory']  + '/_lattice_' + save_file + '.png', dpi=200)
			plt.close()
		return

	def save_numerical_results(self, save_file, eigenenergies):
		with open(self.__configuration['saving_directory']  + '/' + save_file + '.txt', "w") as file:
			for eigenenergy in eigenenergies:
				file.write(str(eigenenergy) +  "\n")
		return

	def get_file_name(self, num_of_atoms):
		sigma = self.__configuration['sigma']
		gauss_sigma = self.__configuration['gauss_sigma']
		x_num_of_steps = self.__configuration['x_num_of_steps']
		calculation_type = self.__configuration['calculation_type']

		if x_num_of_steps != None:
			x_num_of_steps = self.__configuration['x_num_of_steps']
			save_file = calculation_type + \
						'_sigma=' + \
						str(sigma) + \
						'_gauss_sigma=' + \
			 			 str(gauss_sigma) + \
						'_num_of_atoms=' + \
						str(num_of_atoms)
		elif x_num_of_steps == None:
			vertical_num_of_steps = self.__configuration['vertical_num_of_steps']
			horizontal_num_of_steps = self.__configuration['horizontal_num_of_steps']

			save_file = calculation_type + \
			 			'_sigma=' + \
			 			 str(sigma) + \
			 			 '_gauss_sigma=' + \
			 			 str(gauss_sigma) + \
			 			 '_width=' + \
			 			  str(vertical_num_of_steps) + \
			 			  '_length=' + \
			 			  str(horizontal_num_of_steps) + \
			 			  '_num_of_atoms=' + \
			 			  str(num_of_atoms)
		return save_file

