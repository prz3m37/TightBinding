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
        self.__configuration = configuration
        self.__slater_koster = SlaterKoster()
        self.__initialize_tight_binding_params()

    def __initialize_tight_binding_params(self):
        self.power = self.__configuration['power']
        self.epsilon = self.__configuration['epsilon']
        self.distance = self.__configuration['distance']
        self.magnitude = self.__configuration['magnitude']
        self.fermi_level = self.__configuration['fermi_level']
        self.spin_included = self.__configuration['spin_included']
        self.lanczos_vectors = self.__configuration['lanczos_vectors']
        self.number_of_processes = self.__configuration['number_of_cpus']
        self.num_of_eigenvalues = self.__configuration['num_of_eigenvalues']

        self.__set_dimension(spin_included=self.spin_included)

        return

    @staticmethod
    def __set_dimension(spin_included):
        if spin_included == False:
            self.dimension = 10
        else:
            self.dimension = 20

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
        close_friends_indices = ctree.query_ball_tree(
            ctree, r=self.distance, p=self.power, eps=self.epsilon)
        close_friends_indices = list(close_friends_indices)
        return close_friends_indices

    def __calculate_diagonal_elements(self, data, process):
        columns, rows, values = [], [], []
        hosts = data.index
        for host in hosts:
            type_of_host = data['type_of_atom'].iloc[host]
            h_diagonal = self.__slater_koster.calculate_spin_mixing_diagonal()
            h_diagonal_non_zero_indices = np.where(h_diagonal != 0)
            h_diagonal_non_zero_values = np.array(
                h_diagonal[h_diagonal_non_zero_indices])[0]

            h_diagonal_rows = h_diagonal_non_zero_indices[0] + int(
                host * self.dimension)
            h_diagonal_columns = h_diagonal_non_zero_indices[1] + int(
                host * self.dimension)

            columns.append(h_diagonal_columns)
            rows.append(h_diagonal_rows)
            values.append(h_diagonal_non_zero_values)
        return_dict[process] = [columns, rows, values]
        return

    def __calculate_non_diagonal_elements(self, data, close_friends, process):
        columns, rows, values = [], [], []
        for friends in close_friends:
            for friend in friends[1:]:
                host = friends[0]
                type_of_host = data['type_of_atom'].iloc[host]
                ri, rj = data['localization'].iloc[host], data['localization'].iloc[friend]
                type_of_friend = data['type_of_atom'].iloc[friend]
                h_sk = self.__slater_koster.calculate_spin_mixing_sk(
                    type_of_host, type_of_friend)

                h_sk_non_zero_indices = np.where(h_sk != 0)
                h_sk_non_zero_values = np.array(h_sk[h_sk_non_zero_indices])[0]

                h_sk_rows = h_sk_non_zero_indices[0] + \
                    int(host * self.dimension)
                h_sk_columns = h_sk_non_zero_indices[1] + \
                    int(friend * self.dimension)

                columns.append(h_sk_columns)
                rows.append(h_sk_rows)
                values.append(h_sk_non_zero_values)
        return_dict[process] = [columns, rows, values]
        return

    # TODO: Chceck this shit
    def __parallelize_calculations(self, data, close_friends, constants_of_pairs, hosts, atom_store, lp, ld):
        # https://stackoverflow.com/questions/10415028/how-can-i-recover-the-return-value-of-a-function-passed-to-multiprocessing-proce
        number_of_processes = self.__number_of_processes
        manager = multiprocessing.Manager()
        self.__return_dict = manager.dict()
        jobs = []
        for i in range(0, number_of_processes + 1):
            if i == number_of_processes:
                process = multiprocessing.Process(target=self.__calculate_non_diagonal_elements,
                                                  args=(data, hosts, atom_store, lp, ld, i))
            else:
                process = multiprocessing.Process(target=self.__calculate_diagonal_elements,
                                                  args=(data, close_friends, constants_of_pairs, i))
            jobs.append(process)
            process.start()

        for proc in jobs:
            proc.join()
        interaction_matrix_elements = self.return_dict
        return interaction_matrix_elements

    def __run_proper_type_of_calculations(self, **kwargs):
        number_of_processes = self.__number_of_processes
        if number_of_processes == 1 or number_of_processes == None:
            self.__helpers.split_close_friends(
                kwargs, number_of_processes=number_of_processes)
            interaction_matrix_elements = self.__parallelize_calculations(
                kwargs)
        else:
            diag_columns, diag_rows, diag_values = self.__calculate_diagonal_elements(
                kwargs)
            nondiag_columns, nondiag_rows, nondiag_values = self.__calculate_non_diagonal_elements(
                kwargs)
            interaction_matrix_elements = [diag_columns, diag_rows, diag_values,
                                           nondiag_columns, nondiag_rows, nondiag_values]
        return interaction_matrix_elements

    def __construct_interaction_matrix(self, interaction_matrix_elements: (dict, list)):
        if type(interaction_matrix_elements) == dict:
            pass
        if type(interaction_matrix_elements) == list:
            pass
        return

    # TODO split into processes
    def __get_non_zero_values_and_indices(self, data: pd.DataFrame, distance: float, constants_of_pairs: dict, atom_store: dict,
                                          calculation_type: str = 'non spin', method: str = 'distance', number_of_friends: int = None,
                                          lp: float = 0, ld: float = 0, flat: bool = True) -> (list, list, list):
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

        atom_localization = data['localization'].values.tolist()
        close_friends = self.__get_closest_friends(
            atom_localization, method, number_of_friends, distance)

        columns = []
        rows = []
        values = []
        for friends in close_friends:
            host = friends[0]
            type_of_host = data['type_of_atom'].iloc[host]
            h_diagonal = self.__slater_koster.calculate_spin_mixing_diagonal(calculation_type,
                                                                             atom_store,
                                                                             type_of_host,
                                                                             lp,
                                                                             ld)
            h_diagonal_non_zero_indices = np.where(h_diagonal != 0)
            h_diagonal_non_zero_values = np.array(
                h_diagonal[h_diagonal_non_zero_indices])[0]

            h_diagonal_rows = h_diagonal_non_zero_indices[0] + int(
                host * dimension)
            h_diagonal_columns = h_diagonal_non_zero_indices[1] + int(
                host * dimension)

            columns.append(h_diagonal_columns)
            rows.append(h_diagonal_rows)
            values.append(h_diagonal_non_zero_values)
            if len(friends) > 1:
                for friend in friends[1:]:
                    ri, rj = data['localization'].iloc[host], data['localization'].iloc[friend]
                    type_of_friend = data['type_of_atom'].iloc[friend]
                    h_sk = self.__slater_koster.calculate_spin_mixing_sk(calculation_type,
                                                                         ri,
                                                                         rj,
                                                                         constants_of_pairs,
                                                                         type_of_host,
                                                                         type_of_friend,
                                                                         flat)

                    h_sk_non_zero_indices = np.where(h_sk != 0)
                    h_sk_non_zero_values = np.array(
                        h_sk[h_sk_non_zero_indices])[0]

                    h_sk_rows = h_sk_non_zero_indices[0] + \
                        int(host * dimension)
                    h_sk_columns = h_sk_non_zero_indices[1] + \
                        int(friend * dimension)

                    columns.append(h_sk_columns)
                    rows.append(h_sk_rows)
                    values.append(h_sk_non_zero_values)
            else:
                pass

        columns = np.concatenate(columns)
        rows = np.concatenate(rows)
        values = np.concatenate(values)
        self.__helpers.save_log('[INFO]: Interaction matrix calculated \n')

        return columns, rows, values, dimension

    def __create_sparse_matrix(self, number_of_atoms, **kwargs) -> csr_matrix:
        """
        Method converts calculated rows, columns, and non-zero values into sparse matrix.
        Args:
            number_of_atoms: number of atoms
            **kwargs: arguments of __get_non_zero_values_and_indices method
        Returns: sparse interaction matrix for all interacting atoms in lattice in csr format
        """

        columns, rows, values, dimension = self.__get_non_zero_values_and_indices(
            **kwargs)
        matrix_final_1 = coo_matrix((values, (rows, columns)),
                                    shape=(self.dimension * number_of_atoms, self.dimension * number_of_atoms))
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
