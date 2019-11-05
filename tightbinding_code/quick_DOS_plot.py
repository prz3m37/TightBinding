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


  plt.plot(E, D_E)
  plt.axhline(y=0, color='r', linestyle='-')
  plt.xlabel('Energy')
  plt.ylabel('DOS')
  plt.show()

  return

def main():

  start = -10
  stop = 10
  step = 0.01
  gauss_width = 0.07
  file = "/home/przemek/Documents/Modeling/tight_binding/results/non spin_sigma=0.0001_gauss_sigma=0.015_num_of_atoms=4040rectangular.txt"

  quick_plot(file, gauss_width, start, stop, step)

  return

if __name__ == '__main__':
  exit(main())