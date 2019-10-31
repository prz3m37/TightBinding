import numpy as np
import pandas as pd


class LatticeConstructor(object):

    """
    Class responsible for creation of lattice of atoms.
    """

    @staticmethod
    def make_deformation(deforming_func, direction, start, stop):

        size = len(direction)
        deformation = [eval(deforming_func) for x in np.arange(start, stop, (stop - start) / size)]

        for num, loc in enumerate(direction):
            loc[2] = loc[2] + deformation[num]

        return direction

    def construct_skeleton(self, a, vertical_num_of_steps, horizontal_num_of_steps,
                           deformation_func=None, direction=None, start=None, stop=None):

        mother_first_row = [[float(format((-a * np.sqrt(3)/2), '.12g')),float(format( a / 2, '.12g')), 0],
                            [0, float(format( a, '.12g')), 0], 
                            [float(format((a * np.sqrt(3)/2), '.12g')),float(format( a / 2, '.12g')), 0]]

        mother_second_row = [[float(format((-a * np.sqrt(3)/2), '.12g')),float(format( -a / 2, '.12g')), 0],
                             [0, float(format( -a, '.12g')), 0],
                             [float(format((a * np.sqrt(3)/2), '.12g')),float(format( -a / 2, '.12g')), 0]]

        mother_first_row = mother_first_row
        mother_second_row = mother_second_row

        horizontal_shift = np.array([float(format(a * np.sqrt(3), '12g')), 0, 0])
        vertical_shift = np.array([0, float(format(3 * a, '12g')), 0])

        vertical_steps = int(vertical_num_of_steps)
        horizontal_steps = int(horizontal_num_of_steps)

        horizontal_part_first_row = [mother_first_row +(horizontal_step * horizontal_shift)
                                     for horizontal_step in range(0, horizontal_steps)]

        horizontal_part_second_row = [mother_second_row + (horizontal_step * horizontal_shift)
                                      for horizontal_step in range(0, horizontal_steps)]

        horizontal_part_first_row = np.concatenate(horizontal_part_first_row)
        to_drop_from_first_row = []
        for i in range(1, len(horizontal_part_first_row)):
            if np.all(np.isclose(horizontal_part_first_row[i],  horizontal_part_first_row[i-1])) == True:
                to_drop_from_first_row.append(i)
        horizontal_part_first_row_unq = np.delete(horizontal_part_first_row, to_drop_from_first_row, 0)


        horizontal_part_second_row = np.concatenate(horizontal_part_second_row)
        to_drop_from_second_row = []
        for i in range(1, len(horizontal_part_second_row)):
            if np.all(np.isclose(horizontal_part_second_row[i],  horizontal_part_second_row[i-1])) == True:
                to_drop_from_second_row.append(i)
        horizontal_part_second_row_unq = np.delete(horizontal_part_second_row, to_drop_from_second_row, 0)

        horizontal_part_first_row_unq = [list(array) for array in horizontal_part_first_row_unq]
        horizontal_part_second_row_unq = [list(array) for array in horizontal_part_second_row_unq]

        vertical_part_first_row = None
        vertical_part_second_row = None

        if deformation_func is None:

            vertical_part_first_row = [horizontal_part_first_row_unq + (vertical_step * vertical_shift)
                                         for vertical_step in range(0, vertical_steps)]

            vertical_part_second_row = [horizontal_part_second_row_unq + (vertical_step * vertical_shift)
                                          for vertical_step in range(0, vertical_steps)]

        else:

            if direction == 'vertical':
                def_horizontal_first_row = self.make_deformation(deformation_func,
                                                                 horizontal_part_first_row_unq,
                                                                 start,
                                                                 stop)
                def_horizontal_second_row = self.make_deformation(deformation_func,
                                                                  horizontal_part_second_row_unq,
                                                                  start,
                                                                  stop)

                vertical_part_first_row = [def_horizontal_first_row + (vertical_step * vertical_shift)
                                             for vertical_step in range(0, vertical_steps)]

                vertical_part_second_row = [def_horizontal_second_row + (vertical_step * vertical_shift)
                                              for vertical_step in range(0, vertical_steps)]

            elif direction == 'horizontal':
                v_part_first_row = [horizontal_part_first_row_unq + (vertical_step * vertical_shift)
                                    for vertical_step in range(0, vertical_steps)]

                v_part_second_row = [horizontal_part_second_row_unq + (vertical_step * vertical_shift)
                                     for vertical_step in range(0, vertical_steps)]

                vertical_part_first_row = self.make_deformation(deformation_func, v_part_first_row, start, stop)
                vertical_part_second_row = self.make_deformation(deformation_func, v_part_second_row, start, stop)

        lattice = np.concatenate(vertical_part_first_row + vertical_part_second_row)

        return lattice

    @staticmethod
    def add_defects(lattice, number_of_defects):

        random_atoms_to_delete = np.random.randint(low=1, high=lattice.shape[0], size=number_of_defects)
        lattice = np.delete(lattice, random_atoms_to_delete, 0)

        return lattice


    @staticmethod
    def dump_to_dframe(lattice, type_of_atom = 'C'):

        lattice = list(lattice)
        size = len(lattice)
        final_lattice = pd.DataFrame({'number_of_atom': np.arange(0, size, 1), 'localization': lattice,
                                      'type_of_atom': [type_of_atom] * size})
        return final_lattice


