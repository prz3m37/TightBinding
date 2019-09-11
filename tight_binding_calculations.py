import numpy as np
import pandas as pd
from slater_coster_matrix import SlaterKoster
from lattice_maker import LatticeConstructor
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
import sklearn.neighbors as sn
import time
import matplotlib.pyplot as plt


class TightBinding(object):

    """
    Class calculates the matrix of electron interactions between atoms (electronic band structure) in 2D structures
        (for example graphene). Interactions are described in regime of the Tight Binding theorem. Each interaction in
        lattice of atoms is described by Slater-Koster matrix of shape = (10, 10). We take into consideration the
        following orbitals: s, px, py, pz, dxy, dyz, dz2, dxz, dx2dy2, s* in our calculations. The main advantage of
        this program is possibility to calculate electronic band structure for cases like:

        1. Any shape of atom lattices
        2. Lattices with or without defects (i.e. missing atoms)
        3. Multi-level lattices.
        4. Lattices containing different types of atoms.

    """

    def __init__(self) -> None:

        """
        Method calls SlaterKoster and  LatticeConstructor class

        """
        self.__slater_coster = SlaterKoster()
        self.__construct = LatticeConstructor()

    @staticmethod
    def __get_closest_friends(points: list, method: str, number_of_friends: int = None, distance: float = None) -> list:

        """

        Method uses KDTree algorithm for searching closets neighbours (according to Euclidean space) of each atom in
            lattice no matter what shape is.

        Args:
            points: list of coordinates [x,y,z] of all points in lattice.
            method: type of method of searching closest neighbours; if method == 'distance', algorithm will search
                    points in defined distance; else algorithm, will find defined by user number of neighbours.

            number_of_friends: defined by user, number of neighbours of point in lattice
            distance: maximum distance defining point as closest neighbour.

        Returns: List of indices of points which are neighbours.

        """

        tree = sn.KDTree(points, leaf_size=2)
        if method == 'distance':
            close_friends = tree.query_radius(points, r=distance, sort_results=True, return_distance=True)[0]
        else:
            close_friends = tree.query(points, k=number_of_friends)[1]
        return close_friends

    def __get_non_zero_values_and_indices(self, data, calculation_type, method, distance, constants_of_pairs,
                                          atom_store, number_of_friends = None, lp=None, ld=None,
                                          flat=True) -> (list, list, list):

        """

        Method calculates non-zero values (Slater Koster matrix values) with their localization (indices of rows and
            columns) in the final interaction matrix.

        Args:
            data: Input DataFrame where: Indices - integers -standard pd.DataFrame enumeration,
                  Columns: number_of_atom - number of atom in lattice, localization - localization of atom in
                  lattice in format of nested lists [[x1,y1,z1], [x2,y2,z2],...], type_of_atom - type of element of atom
                  string format, for example 'C'

            calculation_type: if calculation_type == 'non spin' Slater Koster matrix will not include Spin-Orbit
                              interactions; else Spin-Orbit interactions will be included

            method: type of method of searching closest neighbours; if method == 'distance', algorithm will search
                    points in defined distance; else algorithm, will find defined by user number of neighbours.

            distance: maximum distance defining point as closest neighbour.
            constants_of_pairs: dict of Slater Koster constants describing interactions between two elements for each
                                orbital; example of input - {('C', 'C'): {'V_sssigma': 0, 'V_spsigma': 0,
                                'V_sdsigma': 0,...

            atom_store: dict of each orbital energies for element; example of input - {'C': {'Es': -8.71, 'Epx': 0,
                        'Epy': 0, ... ; if calculation_type == 'non spin' we have to remember about different spin
                         energies.

            number_of_friends: defined by user, number of neighbours of point in lattice
            lp: p-orbital interaction constant in Spin-Orbit interactions, None if calculation_type == 'non spin'
            ld: d-orbital interaction constant in Spin-Orbit interactions, None if calculation_type == 'non spin'
            flat:

        Returns: Lists of row, column indices and list of non-zero values

        """

        if calculation_type == 'non spin':
            num = 10

        else:
            num = 20

        atom_localization = data['localization'].values.tolist()
        close_friends = self.__get_closest_friends(atom_localization, method, number_of_friends, distance)

        columns = []
        rows = []
        values = []
        for friends in close_friends:

            host = friends[0]
            type_of_host = data['type_of_atom'].iloc[host]
            h_diagonal = self.__slater_coster.calculate_spin_mixing_diagonal(calculation_type,
                                                                             atom_store,
                                                                             type_of_host,
                                                                             lp,
                                                                             ld)
            h_diagonal_non_zero_indices = np.where(h_diagonal != 0)
            h_diagonal_non_zero_values = np.array(h_diagonal[h_diagonal_non_zero_indices])[0]

            h_diagonal_rows = h_diagonal_non_zero_indices[0] + int(host * num)
            h_diagonal_columns = h_diagonal_non_zero_indices[1] + int(host * num)
            
            columns.append(h_diagonal_columns)
            rows.append(h_diagonal_rows)         
            values.append(h_diagonal_non_zero_values)

            if len(friends) != 0:

                for friend in friends[1:]:

                    ri, rj = data['localization'].iloc[host], data['localization'].iloc[friend]
                    type_of_friend = data['type_of_atom'].iloc[friend]

                    h_sk = self.__slater_coster.calculate_spin_mixing_sk(calculation_type,
                                                                         ri,
                                                                         rj,
                                                                         constants_of_pairs,
                                                                         type_of_host,
                                                                         type_of_friend,
                                                                         flat)
                                                                         
                    h_sk_non_zero_indices = np.where(h_sk != 0)
                    h_sk_non_zero_values = np.array(h_sk[h_sk_non_zero_indices])[0]
                    
                    h_sk_rows = h_sk_non_zero_indices[0] + int(host * num)
                    h_sk_columns = h_sk_non_zero_indices[1] + int(friend * num)   
                    
                    columns.append(h_sk_columns)
                    rows.append(h_sk_rows)
                    values.append(h_sk_non_zero_values)

                pass
        
        columns = np.concatenate(columns)
        rows = np.concatenate(rows)
        values = np.concatenate(values)
        print('Matrix calculated')
        
        return columns, rows, values, num

    def __create_sparse_matrix(self, dimension, **kwargs) -> csr_matrix:

        """
        Method converts calculated rows, columns, and non-zero values into sparse matrix.

        Args:
            dimension: dimension of Slater Koster matrix
            **kwargs: arguments of __get_non_zero_values_and_indices method

        Returns: sparse interaction matrix for all interacting atoms in lattice in csr format

        """

        columns, rows, values, num = self.__get_non_zero_values_and_indices(**kwargs)
        matrix_final = coo_matrix((values, (rows, columns)), shape=(dimension * num, dimension * num))
        matrix_final = csr_matrix(matrix_final)

        return matrix_final

    def construct_lattice(self, **kwargs) -> (pd.DataFrame, int):

        """
        Method constructs lattice of atoms and saves it in form of pd.DataFrame
        Args:
            **kwargs: arguments of construct_skeleton method

        Returns: pd.DataFrame with lattice data, and it's shape

        """

        lattice = self.__construct.construct_skeleton(**kwargs)
        lattice_data_frame_format = self.__construct.dump_to_dframe(lattice)
        shape = int(lattice.shape[0])

        return lattice, lattice_data_frame_format, shape

    def calculate_eigenvalues_ang_eigenvectors(self, num_of_eigenvalues, which='LM',
                                               sigma:float=None, **kwargs) -> (np.array, np.array):

        """
        Method calculates eigenvalues (energies) and eigenvectors (wave functions) of interaction matrix.
        Args:
            num_of_eigenvalues: number of eigenvalues to count
            which:
            **kwargs: arguments of __create_sparse_matrix method
            sigma

        Returns: eigenvalues (energies) and eigenvectors (wave functions) of interaction matrix

        """

        sparse_matrix = self.__create_sparse_matrix(**kwargs)
        print('Eigenvalues calculating')
        eigenvalues, eigenvectors = eigsh(sparse_matrix, k=num_of_eigenvalues, which=which, sigma=sigma, mode='normal')

        return eigenvalues, eigenvectors

    @staticmethod
    def evaluate_density_of_states(En, E):

        a = (np.max(En) - np.min(En)) / len(En)
        D_E = [np.sum(np.exp(-(En - e)**2 / (2*a**2)) / (np.pi * a * np.sqrt(2))) for e in E]

        return D_E


atom_store_example1 = {'C': {'Es': -8.7,
                             'Epx': 0,
                             'Epy': 0,
                             'Epz': 0,
                             'Edz2': 0,
                             'Edxz': 0,
                             'Edyz': 0,
                             'Edxy': 0,
                             'Edx2y2': 0,
                             'Estar': 0}}

atom_store_example3 = {'C': {'Es up up': 1,
                             'Epx up up': 1.,
                             'Epy up up': 1,
                             'Epz up up': 1,
                             'Edz2 up up': 1,
                             'Edxz up up': 1,
                             'Edyz up up': 1,
                             'Edxy up up': 1,
                             'Edx2y2 up up': 1,
                             'Estar up up': 1,
                             'Es down down': 1,
                             'Epx down down': 1,
                             'Epy down down': 1,
                             'Epz down down': 1,
                             'Edz2 down down': 1,
                             'Edxz down down': 1,
                             'Edyz down down': 1,
                             'Edxy down down': 1,
                             'Edx2y2 down down': 1,
                             'Estar down down': 1}}


mapping_system = pd.DataFrame({'number_of_atom': [0, 1, 2, 3, 4, 5],
                               'type_of_atom': ['C', 'C', 'C', 'C', 'C', 'C'],
                               'localization': [np.array([-np.sqrt(3)/2, -0.5, 0]),
                                                np.array([0, -1, 0]),
                                                np.array([np.sqrt(3)/2, -0.5, 0]),
                                                np.array([-np.sqrt(3)/2, 0.5, 0]),
                                                np.array([0, 1, 0]),
                                                np.array([np.sqrt(3)/2, 0.5, 0])]})


constans_of_pairs_example1 = {('C', 'C'): {'V_sssigma': -6.7,
                                           'V_spsigma': 5.5,
                                           'V_sdsigma': 0,
                                           'V_starssigma': 0,
                                           'V_starpsigma': 0,
                                           'V_stardsigma': 0,
                                           'V_ppsigma': 5.1,
                                           'V_pppi': -3.1,
                                           'V_pdsigma': 0,
                                           'V_pdpi': 0,
                                           'V_ddsigma': 0,
                                           'V_ddpi': 0,
                                           'V_ddd': 0}}

# TODO taśma grafenowa idealna , a potem robie dziury, parametryzacja dla grafenu zapytac google tight binding. gdy nie mam atomu potencjał bardzo duzy sprawdzenie on site energies bardzo duze -  wygenerowac
# TODO 5000 tys atomów 


start = time.time()
tb = TightBinding()
lc = LatticeConstructor()
lattice, lattice_data_frame_format, dims = tb.construct_lattice(a=1.42,
                                                                vertical_num_of_atoms=50,
                                                                horizontal_num_of_atoms=5)
lc.plot_lattice(lattice)

print('Calculation for', str(dims), ' atoms')
sigma = 0.5
energy, wave_function = tb.calculate_eigenvalues_ang_eigenvectors(dimension=dims,
                                                                  data=lattice_data_frame_format,
                                                                  calculation_type='non spin',
                                                                  method='distance',
                                                                  distance=1,
                                                                  constants_of_pairs=constans_of_pairs_example1,
                                                                  atom_store=atom_store_example1,
                                                                  num_of_eigenvalues=60,
                                                                  which='BE',
                                                                  sigma=sigma)
E = np.arange(-20, 20, 0.1)

# energy = sigma + (1 / energy)
density_of_states = tb.evaluate_density_of_states(energy, E)


end = time.time()
print('Calculation time for ' + str(dims) + ' atoms: ', (end - start) / 60., ' minutes')


print(energy)
plt.plot(E, density_of_states)
plt.show()
