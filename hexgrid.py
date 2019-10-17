import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt


def generate_middle_strip(a, x_num_of_steps):

    mother_first_row = [[-a * np.sqrt(3)/2, a / 2, 0], [0, a, 0], [a * np.sqrt(3)/2, a / 2, 0]]
    mother_second_row = [[-a * np.sqrt(3)/2, -a / 2, 0], [0, -a, 0], [a * np.sqrt(3)/2, -a / 2, 0]]

    x_step = np.array([a * np.sqrt(3), 0, 0])
    x_num_of_steps = int(x_num_of_steps)

    x_part_first_row = [mother_first_row + (x_num * x_step) for x_num in range(0, x_num_of_steps)]
    x_part_second_row = [mother_first_row + (x_num * x_step) for x_num in range(0, x_num_of_steps)]

    x_part_first_row = np.concatenate(x_part_first_row)
    to_drop_from_first_row = []
    for i in range(1, len(x_part_first_row)):
        if np.all(np.isclose(x_part_first_row[i],  x_part_first_row[i-1])) == True:
            to_drop_from_first_row.append(i)
    x_part_first_row_unq = np.delete(x_part_first_row, to_drop_from_first_row, 0)
        

    x_part_second_row = np.concatenate(x_part_second_row)
    to_drop_from_second_row = []
    for i in range(1, len(x_part_second_row)):
        if np.all(np.isclose(x_part_second_row[i],  x_part_second_row[i-1])) == True:
            to_drop_from_second_row.append(i)
    x_part_second_row_unq = np.delete(x_part_second_row, to_drop_from_second_row, 0)

    x_part_first_row_unq = [list(array) for array in x_part_first_row_unq]
    x_part_second_row_unq = [list(array) for array in x_part_second_row_unq]

    '''
    x_part_first_row_unq = np.concatenate(x_part_first_row_unq)
    x_part_second_row_unq = np.concatenate(x_part_second_row_unq)
    '''
    return x_part_first_row_unq, x_part_second_row_unq
    
    
def generate_hexgrid(a, x_num_of_steps):
    
    y_step = np.array([0, 3 * a, 0])
    y_num_of_steps = 0
    for num in range(3, x_num_of_steps + 1):  #CHYBA OK ALE POPATRZEC JESZCZE
        if num % 2 != 0:
            y_num_of_steps += 1 
    
    first_row, second_row = generate_middle_strip(a, x_num_of_steps)
    
    y_part_first_row = [] 
    for num, y_num in enumerate(range(0, y_num_of_steps)):
        next_row = first_row + (y_num * y_step)
        y_part_first_row.append(next_row[num:-num])
        
    
    y_part_second_row = [] 
    for num, y_num in enumerate(range(0, y_num_of_steps)):
        next_row = second_row + (y_num * y_step)
        y_part_second_row.append(next_row[num:-num])
    
    hexgrid = np.concatenate(y_part_first_row + y_part_second_row + first_row + second_row)
    return hexgrid
    
    
def dump_to_dframe(lattice):

    lattice = list(lattice)
    size = len(lattice)
    final_lattice = pd.DataFrame({'number_of_atom': np.arange(0, size, 1), 'localization': lattice,
                                  'type_of_atom': ['C'] * size})
    return final_lattice
    

def plot_lattice(lattice, projection=None):

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
        plt.show()
    else:

        plt.scatter(x, y, marker='o')
        plt.show()

    return
