import math
import decimal
import numpy as np
import pandas as pd 

class HexGrid():

    def generate_middle_strip(self, a, x_num_of_steps):
        '''
        num = decimal.Decimal(str(a))
        decimals = abs(num.as_tuple().exponent)
        decimals = 12
        mother_first_row = [[float(format((-a * np.sqrt(3)/2), '.%sg' % (decimals))),float(format( a / 2, '.%sg' % (decimals))), 0],
                            [0, float(format( a, '.%sg' % (decimals))), 0], 
                            [float(format((a * np.sqrt(3)/2), '.%sg' % (decimals))),float(format( a / 2, '.%sg' % (decimals))), 0]]

        mother_second_row = [[float(format((-a * np.sqrt(3)/2), '.%sg' % (decimals))),float(format( -a / 2, '.%sg' % (decimals))), 0],
                             [0, float(format( -a, '.%sg' % (decimals))), 0],
                             [float(format((a * np.sqrt(3)/2), '.%sg' % (decimals))),float(format( -a / 2, '.%sg' % (decimals))), 0]]
        '''

        mother_first_row = [[-a * np.sqrt(3)/2, a / 2, 0],
                            [0, a, 0], 
                            [a * np.sqrt(3)/2, a / 2, 0]]

        mother_second_row = [[-a * np.sqrt(3)/2, -a / 2, 0],
                            [0, -a, 0], 
                            [a * np.sqrt(3)/2, -a / 2, 0]]

        x_step = np.array([a * np.sqrt(3), 0, 0])
        x_num_of_steps = int(x_num_of_steps)

        x_part_first_row = [mother_first_row + (x_num * x_step) for x_num in range(0, x_num_of_steps)]
        x_part_second_row = [mother_second_row + (x_num * x_step) for x_num in range(0, x_num_of_steps)]

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
        
        
    def generate_hexgrid(self, a, x_num_of_steps):
        
        first_row, second_row = self.generate_middle_strip(a, x_num_of_steps)
        y_step = np.array([0, 3 * a, 0])
        y_num_of_steps = 0
        for num in range(3, x_num_of_steps + 1):
            if num % 2 != 0:
                y_num_of_steps += 1
        
        #number_of_atoms = np.sum([(2 * x_num_of_steps + 1) - (2 * step) for step in range(0, y_num_of_steps+1)]) 
        #print('number_of_atoms', 2 * number_of_atoms)
        
        if y_num_of_steps % 2 == 0:
           y_red = int(y_num_of_steps / 2)
           y_blue = y_red
        else:
            y_red = int(y_num_of_steps / 2) + 1
            y_blue = y_red - 1
            
        y_up_part = []
        y_num_red_temp = 1
        for y_num_red in range(1, y_red + 1):
            next_first_row = second_row + (y_num_red * y_step)
            y_up_part.append(next_first_row[y_num_red_temp:-y_num_red_temp])
            y_num_red_temp = y_num_red_temp + 2

        for y_num_blue in range(1, y_blue + 1):
            next_second_row = first_row + (y_num_blue * y_step)
            y_up_part.append(next_second_row[y_num_blue*2:-y_num_blue*2])
            
        y_down_part = [] 
        for y_num_blue in range(1, y_blue + 1):
            next_first_row = second_row - (y_num_blue * y_step)
            y_down_part.append(next_first_row[2*y_num_blue:-2*y_num_blue])
              
        y_num_red_temp = 1
        for y_num_red in range(1, y_red + 1):
            next_second_row = first_row - (y_num_red * y_step)
            y_down_part.append(next_second_row[y_num_red_temp:-y_num_red_temp])
            y_num_red_temp = y_num_red_temp + 2
        
        y_part_first_row= np.concatenate(y_up_part)
        y_part_second_row= np.concatenate(y_down_part)
        
        y_part_first_row = [array.tolist() for array in y_part_first_row]
        y_part_second_row = [array.tolist() for array in y_part_second_row]
               
        hexgrid = y_part_first_row + y_part_second_row + first_row + second_row
        hexgrid = [np.array(l) for l in hexgrid]
        return hexgrid
        
    @staticmethod    
    def dump_to_dframe(lattice):

        lattice = list(lattice)
        size = len(lattice)
        final_lattice = pd.DataFrame({'number_of_atom': np.arange(0, size, 1), 'localization': lattice,
                                      'type_of_atom': ['C'] * size})

        #remove_n = int(size * 0.4)
        #drop_indices = np.random.choice(final_lattice.index, remove_n, replace=False)
        #final_lattice_dropped = final_lattice.drop(drop_indices)
        return final_lattice
        







