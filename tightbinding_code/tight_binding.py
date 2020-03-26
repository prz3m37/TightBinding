import gc
import numpy as np
import pandas as pd
import multiprocessing
import scipy.spatial as scsp
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

    def __init__(self, helpers, configuration) -> None:
        """
        Method calls SlaterKoster class
        """
        self.__return_dict = None
        self.__helpers = helpers
        self.__params = configuration
        self.__load_tight_binding_params()
        self.__slater_koster = SlaterKoster(self.spin_included)

    def __initialize_tight_binding_params(self):
        self.power = self.__params['power']
        self.epsilon = self.__params['epsilon']
        self.distance = self.__params['distance']
        self.magnitude = self.__params['magnitude']
        self.fermi_level = self.__params['fermi_level']
        self.spin_included = self.__params['spin_included']
        self.lanczos_vectors = self.__params['lanczos_vectors']
        self.number_of_processes = self.__params['number_of_cpus']
        self.num_of_eigenvalues = self.__params['num_of_eigenvalues']

        self.__set_dimension(spin_included=self.spin_included)

        return

    @staticmethod
    def __set_dimension(spin_included):
        if spin_included == False:
            self.dimension = 10
        else:
            self.dimension = 20

    @staticmethod
    def __get_non_zero_indices_and_values(dimension, matrix, atom_i, atom_j):

        matrix_non_zero_indices = np.where(matrix != 0)
        matrix_non_zero_values = np.array(matrix[matrix_non_zero_indices])[0]
        matrix_rows = matrix_non_zero_indices[0] + \
            int(atom_i * dimension)
        matrix_columns = matrix_non_zero_indices[1] + \
            int(atom_j * dimension)

        return matrix_columns, matrix_rows, matrix_non_zero_values

    # TODO: Change for scipy version
    def __get_closest_friends(self, points: list) -> list:
        """
        Method uses KDTree algorithm for searching closets neighbours (according to Euclidean space) of each atom in
            lattice no matter what shape is.
        Args:
            points: list of coordinates [x,y,z] of all points in lattice.
            distance: maximum distance defining point as closest neighbour.
            power: Minkowski norm
            epsilon: Approximate search. Branches of the tree are not explored if their nearest points are further
             than r/(1+eps), and branches are added in bulk if their furthest points are nearer than r * (1+eps).
              eps has to be non-negative.
        Returns: List of indices of points which are neighbours.
        """
        ctree = scsp.cKDTree(points)
        close_friends_indices = ctree.query_pairs(
            r=self.distance, p=self.power, eps=self.epsilon)
        close_friends_indices = list(close_friends_indices)
        return close_friends_indices

    def __calculate_diagonal_elements(self, data: pd.DataFrame, process: (int, None) = None):
        columns, rows, values = [], [], []
        hosts = data.index
        for host in hosts:
            type_of_host = data['type_of_atom'].iloc[host]
            h_diagonal = self.__slater_koster.calculate_spin_mixing_diagonal(
                type_of_host)

            h_diagonal_columns, h_diagonal_rows, h_diagonal_non_zero_values = \
                self.__get_non_zero_indices_and_values(
                    matrix=h_diagonal, atom_i=host, atom_j=host, dimension=self.dimension)

            columns.append(h_diagonal_columns)
            rows.append(h_diagonal_rows)
            values.append(h_diagonal_non_zero_values)
        return_dict[process] = [columns, rows, values]
        return

    def __calculate_non_diagonal_elements(self, data: pd.DataFrame, close_friends: np.array(),
                                          process: (int, None) = None):
        columns, rows, values = [], [], []
        for friends in close_friends:
            for friend in friends[1:]:
                host = friends[0]
                type_of_host = data['type_of_atom'].iloc[host]
                ri, rj = data['localization'].iloc[host], data['localization'].iloc[friend]
                type_of_friend = data['type_of_atom'].iloc[friend]
                h_sk = self.__slater_koster.calculate_spin_mixing_sk(
                    ri=ri, rj=rj, atom_i=type_of_host, atom_j=type_of_friend)

            h_sk_columns_ij, h_sk_rows_ij, h_sk_non_zero_values_ij = \
                self.__get_non_zero_indices_and_values(
                    matrix=h_sk, atom_i=host, atom_j=friend, dimension=self.dimension)
            h_sk_columns_ji, h_sk_rows_ji, h_sk_non_zero_values_ji = \
                self.__get_non_zero_indices_and_values(
                    matrix=h_sk, atom_i=friend, atom_j=host, dimension=self.dimension)

            h_sk_columns = np.concatenate((h_sk_columns_ij, h_sk_columns_ji))
            h_sk_rows = np.concatenate((h_sk_rows_ji, h_sk_rows_ji))
            h_sk_non_zero_values = np.concatenate(
                (h_sk_non_zero_values_ij, h_sk_non_zero_values_ji))

            columns.append(h_sk_columns)
            rows.append(h_sk_rows)
            values.append(h_sk_non_zero_values)

        return_dict[process] = [columns, rows, values]
        return

    def __parallelize_calculations(self, data: pd.DataFrame, splitted_close_friends: np.array,
                                   number_of_processes: (int, None)):
        manager = multiprocessing.Manager()
        self.__return_dict = manager.dict()
        jobs = []
        for i, close_friends in zip(range(0, number_of_processes + 1), splitted_close_friends):
            if i == number_of_processes:
                process = multiprocessing.Process(target=self.__calculate_non_diagonal_elements,
                                                  args=(data, i))
            else:
                process = multiprocessing.Process(target=self.__calculate_diagonal_elements,
                                                  args=(data, close_friends, i))
            jobs.append(process)
            process.start()

        for proc in jobs:
            proc.join()
        interaction_matrix_elements = self.return_dict
        return interaction_matrix_elements

    def __run_interaction_matrix_calculations(self, data: pd.DataFrame):
        atom_localizations = data['localization'].values.tolist()
        close_friends = self.__get_closest_friends(points=atom_localizations)
        if self.__number_of_processes == 1 or self.__number_of_processes == None:
            number_of_processes = 1
            diag_columns, diag_rows, diag_values = \
                self.__calculate_diagonal_elements(data=data)
            nondiag_columns, nondiag_rows, nondiag_values = \
                self.__calculate_non_diagonal_elements(
                    data=data, close_friends=close_friends)
            interaction_matrix_elements = [diag_columns, nondiag_columns,
                                           diag_rows, nondiag_rows,
                                           diag_values, nondiag_values]
        else:
            splitted_close_firends = self.__helpers.split_close_friends(
                close_friends=close_friends, number_of_cpus=number_of_processes)
            interaction_matrix_elements = \
                self.__parallelize_calculations(
                    data=data, splitted_close_friends=splitted_close_firends,
                    number_of_processes=number_of_processes)

        self.__helpers.save_log(
            '[INFO]: Calculations are runned on {} processes/CPUS\n'.format(number_of_processes))
        return interaction_matrix_elements

    def __construct_interaction_matrix(self, data: pd.DataFrame):
        interaction_matrix_elements = self.__run_interaction_matrix_calculations(
            data)
        if self.__number_of_processes == 1 or self.__number_of_processes == None:
            columns = np.concatenate(interaction_matrix_elements[0:2])
            rows = np.concatenate(interaction_matrix_elements[2:4])
            values = np.concatenate(interaction_matrix_elements[4:6])
        else:
            columns_contener, rows_contener, values_contener = [], [], []
            for i in range(0, number_of_processes + 1):
                columns, rows, values = interaction_matrix_elements[i]
                columns_contener.append(columns)
                rows_contener.append(rows)
                values_contener.append(values)

            columns = np.concatenate(columns_contener)
            rows = np.concatenate(rows_contener)
            values = np.concatenate(values_contener)

        self.__helpers.save_log('[INFO]: Interaction matrix calculated \n')
        return columns, rows, values

    def __create_sparse_matrix(self, number_of_atoms: int, **kwargs) -> csr_matrix:
        """
        Method converts calculated rows, columns, and non-zero values into sparse matrix.
        Args:
            number_of_atoms: number of atoms
            **kwargs: arguments of __get_non_zero_values_and_indices method
        Returns: sparse interaction matrix for all interacting atoms in lattice in csr format
        """

        columns, rows, values = self.__construct_interaction_matrix(
            **kwargs)
        matrix_final_1 = coo_matrix((values, (rows, columns)),
                                    shape=(self.dimension * number_of_atoms,
                                           self.dimension * number_of_atoms))
        matrix_final = csr_matrix(matrix_final_1)

        return matrix_final

    def calculate_eigenvalues_and_eigenvectors(self, **kwargs) -> (np.array, np.array):
        """
        Method calculates eigenvalues (energies) and eigenvectors (wave functions) of interaction matrix.
        Args:
            num_of_eigenvalues: number of eigenvalues to count
            magnitude: type of evalueted eigenvalues.
            **kwargs: arguments of __create_sparse_matrix method
            fermi_level: shift-invert parameter (in fact this can be interpreted as Fermi energy level)
        Returns: eigenvalues (energies) and eigenvectors (wave functions) of interaction matrix
        """

        interaction_matrix = self.__create_sparse_matrix(**kwargs)
        self.__helpers.save_log(
            '[INFO]: Eigenvalues calculating from matrix of shape: ', interaction_matrix.shape, '\n')
        eigenvalues, eigenvectors = eigsh(interaction_matrix,
                                          k=self.num_of_eigenvalues,
                                          which=self.magnitude,
                                          sigma=self.fermi_level,
                                          mode='normal',
                                          ncv=self.lanczos_vectors)

        return eigenvalues, eigenvectors, interaction_matrix

    def evaluate_density_of_states(self, eigenenergies: np.array, E: np.array, gauss_width: float) -> list:
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
            D_E = D_E + np.exp(-(E - eigenenergy)**2 / (2 *
                                                        gauss_width**2)) / (np.pi * gauss_width * np.sqrt(2))
        self.__helpers.save_log('[INFO]: DOS calculated \n')
        return D_E

    def evaluate_projected_density_of_states(self, eigenenergies: np.array, E: np.array, eigenstates: np.array,
                                             gauss_width: float) -> list:
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
