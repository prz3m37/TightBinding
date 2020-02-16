import gc
import numpy as np
import pandas as pd
import sklearn.neighbors as sn
from scipy.sparse.linalg import eigsh
from slater_koster import SlaterKoster
from scipy.sparse import coo_matrix, csr_matrix


class TightBinding(object):

    """
    Class calculates the matrix of electron interactions between atoms (electronic band structure) in 2D structures
        (for example graphene). Interactions are described in regime of the Tight Binding theorem. Each interaction in
        lattice of atoms is described by Slater-Koster matrix of shape = (10, 10). We take into consideration the
        following orbitals: s, px, py, pz, dxy, dyz, dxz, dx2dy2, dz2, s* in our calculations. The main advantage of
        this program is possibility to calculate electronic band structure for cases like:
        1. Any shape of atom lattices
        2. Lattices with or without defects (i.e. missing atoms)
        3. Multi-level lattices.
        4. Lattices containing different types of atoms.
    """

    def __init__(self, helpers) -> None:

        """
        Method calls SlaterKoster class
        """
        self.__helpers = helpers
        self.__slater_coster = SlaterKoster()

    @staticmethod
    def __get_closest_friends(points: list, method: str, number_of_friends:int=None, distance: float = None) -> list:

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

        tree = sn.KDTree(np.array(points), leaf_size=2)
        if method == 'distance':
            close_friends, _ = tree.query_radius(points, r=distance, sort_results=True, return_distance=True)
        else:
            close_friends = tree.query(points, k=number_of_friends)[1]

        return close_friends

    def __get_non_zero_values_and_indices(self, data: pd.DataFrame, distance: float, constants_of_pairs: dict, atom_store: dict, 
                                          calculation_type: str ='non spin', method: str = 'distance', number_of_friends: int=None, 
                                          lp: float=0, ld: float=0, flat: bool=True) -> (list, list, list):

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
            if len(friends) > 1:
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
            else:
                pass

        columns = np.concatenate(columns)
        rows = np.concatenate(rows)
        values = np.concatenate(values)
        self.__helpers.save_log('[INFO]: Interaction matrix calculated \n')
        
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
        matrix_final_1 = coo_matrix((values, (rows, columns)), shape=(dimension * num, dimension * num))
        matrix_final = csr_matrix(matrix_final_1)

        return matrix_final
        
    def calculate_eigenvalues_and_eigenvectors(self, num_of_eigenvalues: int, magnitude: str = 'LM',
                                               fermi_level: float=None, lanczos_vectors=None, **kwargs) -> (np.array, np.array):

        """
        Method calculates eigenvalues (energies) and eigenvectors (wave functions) of interaction matrix.
        Args:
            num_of_eigenvalues: number of eigenvalues to count
            magnitude: type of evalueted eigenvalues.
            **kwargs: arguments of __create_sparse_matrix method
            fermi_level: shift-invert parameter (in fact this can be interpreted as Fermi energy level)
        Returns: eigenvalues (energies) and eigenvectors (wave functions) of interaction matrix
        """

        sparse_matrix = self.__create_sparse_matrix(**kwargs)
        self.__helpers.save_log('[INFO]: Eigenvalues calculating from matrix of shape: ', sparse_matrix.shape, '\n')
        eigenvalues, eigenvectors = eigsh(sparse_matrix, 
                                          k=num_of_eigenvalues, 
                                          which=magnitude, 
                                          sigma=fermi_level, 
                                          mode='normal', 
                                          ncv=lanczos_vectors)

        return eigenvalues, eigenvectors

    def diagonalize_tb_matrix(self, sparse_matrix, num_of_eigenvalues: int, magnitude: str = 'LM',
                              fermi_level: float=None, lanczos_vectors=None):

        self.__helpers.save_log('[INFO]: Eigenvalues calculating from matrix of shape: ', sparse_matrix.shape, '\n')
        eigenvalues, eigenvectors = eigsh(sparse_matrix,
                                          k=num_of_eigenvalues,
                                          which=magnitude,
                                          sigma=fermi_level,
                                          mode='normal',
                                          ncv=lanczos_vectors)

        return eigenvalues, eigenvectors

    #@staticmethod
    def evaluate_density_of_states(self, eigenenergies:np.array, E:np.array, gauss_width:float) -> list:
        """
        Method calculates density of states. Delta function is approximated by Gauss function.
        Args:
            eigenenergies: eigenvalues of interaction matrix
            E: list of arguments
            gauss_width:
        Returns: Density of states
        """

        D_E = 0
        for eigenenergy in eigenenergies:
            D_E = D_E + np.exp(-(E - eigenenergy)**2 / (2 * gauss_width**2)) / (np.pi * gauss_width * np.sqrt(2))
        self.__helpers.save_log('[INFO]: DOS calculated \n')
        return D_E

    #@staticmethod
    def evaluate_projected_density_of_states(self, eigenenergies:np.array, E:np.array, eigenstates:np.array,
                                             gauss_width:float) -> list:

        """
        Method calculates projected density of states
        Args:
            eigenenergies: eigenvalues of interaction matrix
            E: list of arguments
            eigenstates: eigenvectors of interaction matrix
            gauss_width: Gauss function width
        Returns: Projected density of states
        """

        N = len(eigenenergies)
        a = (np.max(eigenenergies) - np.min(eigenenergies)) / N
        D_projected = 0
        for num, eigenenergy in enumerate(eigenenergies):
            D_projected = D_projected + np.abs(np.dot(np.conj(eigenstates[:, num]),
                                                      eigenstates[:, np.where(eigenenergy)[0][0]]))\
                          * np.exp(-(E - eigenenergy)**2 / (2 * gauss_width**2)) / (np.pi * gauss_width * np.sqrt(2))
        self.__helpers.save_log('[INFO]: Projected DOS calculated \n')
        return D_projected
