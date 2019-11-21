import numpy as np
import matplotlib.pyplot as plt

def quick_plot(results_file, gauss_width, start, stop, step):

  with open(results_file, "r") as results:
    results = results.read().split('\n')
    results = [float(res) for res in results[:-1]]

  eigenenergies = results
  gauss_width = gauss_width

  D_E = 0
  E = np.arange(start, stop, step)
  for eigenenergy in eigenenergies:
      D_E = D_E + np.exp(-(E - eigenenergy)**2 / (2 * gauss_width**2)) / (np.pi * gauss_width * np.sqrt(2))

  plt.figure(figsize=(13.66, 7.68))
  plt.plot(E, D_E)
  plt.axhline(y=0, color='r', linestyle='-')
  plt.xlabel('Energy [a.u.]')
  plt.ylabel('DOS')
  #plt.title('Density of states for ' + str(num_of_atoms) + ' atoms')
  #plt.savefig(self.__configuration['saving_directory'] + '/__DOS__' + save_file + '.png', dpi=400)
  plt.show()
  return

def main():

  start = -18
  stop = 18
  step = 0.01
  gauss_width = 0.04
  path = "/home/przemek/Documents/Modeling/tight_binding/results_diploma/"

  results_file = ['2019-11-19 11:40:31.973504_non spin_sigma=0.0001_gauss_sigma=0.015_num_of_atoms=5100_rectagonal1.txt',
  '2019-11-19 13:46:42.980739_non spin_sigma=0.0001_gauss_sigma=0.015_num_of_atoms=5100_rectagonal2.txt',
  '2019-11-20 01:47:24.023660_non spin_sigma=0.0001_gauss_sigma=0.015_num_of_atoms=4080_rectagonal2.txt',
  '2019-11-20 02:59:03.333819_non spin_sigma=0.0001_gauss_sigma=0.015_num_of_atoms=4080_rectagonal3.txt',
  '2019-11-20 04:07:52.628016_non spin_sigma=0.0001_gauss_sigma=0.015_num_of_atoms=4840_rectagonal4.txt',
  '2019-11-20 05:28:16.494373_non spin_sigma=0.0001_gauss_sigma=0.015_num_of_atoms=4840_rectagonal5.txt',
  '2019-11-20 06:50:22.335714_non spin_sigma=0.0001_gauss_sigma=0.015_num_of_atoms=4840_rectagonal6.txt',
  '2019-11-20 08:40:02.753147_non spin_sigma=0.0001_gauss_sigma=0.015_num_of_atoms=4840_rectagonal7.txt',
  '2019-11-20 11:27:08.082522_non spin_sigma=0.0001_gauss_sigma=0.015_num_of_atoms=4080_rectagonal8.txt']

  for result_file in results_file:
    quick_plot(path + result_file, gauss_width, start, stop, step)

  return

if __name__ == '__main__':
  exit(main())

