import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class LatticeConstructor(object):

    @staticmethod
    def make_deformation(deforming_func, direction, start, stop):

        size = len(direction)
        deformation = [eval(deforming_func) for x in np.arange(start, stop, (stop - start) / size)]

        for num, loc in enumerate(direction):
            loc[2] = loc[2] + deformation[num]

        return direction

    def construct_skeleton(self, a, vertical_num_of_atoms, horizontal_num_of_atoms,
                           deformation_func=None, direction=None, start=None, stop=None):

        mother_first_row = [[-a * np.sqrt(3)/2, a / 2, 0], [0, a, 0], [a * np.sqrt(3)/2, a / 2, 0]]
        mother_second_row = [[-a * np.sqrt(3)/2, -a / 2, 0], [0, -a, 0], [a * np.sqrt(3)/2, -a / 2, 0]]

        vertical_shift = np.array([a * np.sqrt(3), 0, 0])
        horizontal_shift = np.array([0, 3 * a, 0])

        vertical_steps = int(vertical_num_of_atoms)
        horizontal_steps = int(horizontal_num_of_atoms)

        vertical_part_first_row = [mother_first_row + (vertical_step * vertical_shift)
                                   for vertical_step in range(0, vertical_steps)]

        vertical_part_first_row_unq = {array.tostring(): list(array)
                                       for array in np.concatenate(vertical_part_first_row)}

        vertical_part_first_row_unq = list(vertical_part_first_row_unq.values())

        vertical_part_second_row = [mother_second_row + (vertical_step * vertical_shift)
                                    for vertical_step in range(0, vertical_steps)]

        vertical_part_second_row_unq = {array.tostring(): list(array)
                                        for array in np.concatenate(vertical_part_second_row)}

        vertical_part_second_row_unq = list(vertical_part_second_row_unq.values())

        horizontal_part_first_row = None
        horizontal_part_second_row = None

        if deformation_func is None:

            horizontal_part_first_row = [vertical_part_first_row_unq + (horizontal_step * horizontal_shift)
                                         for horizontal_step in range(0, horizontal_steps)]

            horizontal_part_second_row = [vertical_part_second_row_unq + (horizontal_step * horizontal_shift)
                                          for horizontal_step in range(0, horizontal_steps)]

        else:

            if direction == 'vertical':
                def_vertical_first_row = self.make_deformation(deformation_func,
                                                               vertical_part_first_row_unq,
                                                               start,
                                                               stop)
                def_vertical_second_row = self.make_deformation(deformation_func,
                                                                vertical_part_second_row_unq,
                                                                start,
                                                                stop)

                horizontal_part_first_row = [def_vertical_first_row + (horizontal_step * horizontal_shift)
                                             for horizontal_step in range(0, horizontal_steps)]

                horizontal_part_second_row = [def_vertical_second_row + (horizontal_step * horizontal_shift)
                                              for horizontal_step in range(0, horizontal_steps)]

            elif direction == 'horizontal':
                h_part_first_row = [vertical_part_first_row_unq + (horizontal_step * horizontal_shift)
                                    for horizontal_step in range(0, horizontal_steps)]

                h_part_second_row = [vertical_part_second_row_unq + (horizontal_step * horizontal_shift)
                                     for horizontal_step in range(0, horizontal_steps)]

                horizontal_part_first_row = self.make_deformation(deformation_func, h_part_first_row, start, stop)
                horizontal_part_second_row = self.make_deformation(deformation_func, h_part_second_row, start, stop)

        lattice = np.concatenate(horizontal_part_first_row + horizontal_part_second_row)

        return lattice

    @staticmethod
    def dump_to_dframe(lattice):

        lattice = list(lattice)
        size = len(lattice)
        final_lattice = pd.DataFrame({'number_of_atom': np.arange(0, size, 1), 'localization': lattice,
                                      'type_of_atom': ['C'] * size})
        return final_lattice

    @staticmethod
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


