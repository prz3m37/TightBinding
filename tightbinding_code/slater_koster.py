import numpy as np


class SlaterKoster(object):

    """
    SlaterKoster class calculates Slater-Koster matrix which describes interactions between two atoms in lattice.
        User can define whether wants to take Spin Orbit interaction into consideration or not.
    """

    @staticmethod
    def get_atom_type(atom_store: dict, atom_type: str) -> (dict, None):

        """
        Method passes dictionary of physical constants of atoms in lattice
        Args:
            atom_store: dict with physical constants (energies of bands) of atoms in lattice
            atom_type: type of element (for example 'C')
        Returns: dict  dict with physical constants (energies of bands) of atoms in lattice
        """

        if atom_store is None:
            atom = None
        else:
            atom = atom_store[atom_type]

        return atom

    @staticmethod
    def __get_interaction_constants(constants: dict, atom_i_type: str, atom_j_type: str) -> (dict, None):

        """
        Method returns dictionary with interaction constants. If constants are None,
            then all constants are equal to zero
        Args:
            constants: interaction constants
            atom_i_type: type of element interacting atom on i-th position
            atom_j_type: type of element interacting atom on j-th position
        Returns: dictionary with interaction constants or None if constants == None
        """

        if constants is None:

            const = {'V_sssigma': 0,
                     'V_spsigma': 0,
                     'V_sdsigma': 0,
                     'V_starsssigma':0,
                     'V_starssigma': 0,
                     'V_starpsigma': 0,
                     'V_stardsigma': 0,
                     'V_ppsigma': 0,
                     'V_pppi': 0,
                     'V_pdsigma': 0,
                     'V_pdpi': 0,
                     'V_ddsigma': 0,
                     'V_ddpi': 0,
                     'V_ddd': 0}
        else:
            const = constants[atom_i_type, atom_j_type]

        return const

    @staticmethod
    def __get_directional_cosines(ri: np.array, rj: np.array) -> dict:

        """
        Method calculates directional cosines for atoms in lattice
        Args:
            ri: position of i-th atom in lattice
            rj: position of j-th atom in lattice
        Returns: directional cosines
        """
        r_diff = ri - rj
        r_diff_mod = np.linalg.norm(ri - rj)
        n = {'nx': r_diff[0] / r_diff_mod, 'ny': r_diff[1] / r_diff_mod, 'nz': r_diff[2] / r_diff_mod}
        return n

    def __calculate_slayterkoster_matrix(self, ri: int, rj: int,
                                         constants_of_pairs: dict, atom_i_type: str, atom_j_type: str) -> np.matrix:

        """
        Method calculates Slater-Koster matrix (dimension: 10x10)
        Args:
            ri: position of i-th atom in lattice
            rj: position of j-th atom in lattice
            constants_of_pairs: dictionary with interaction constants
            atom_i_type: type of element i-th atom
            atom_j_type: type of element j-th atom
        Returns: Slater-Koster matrix
        """

        H_SK = np.zeros((10, 10))
        n = self.__get_directional_cosines(ri, rj)
        integral_const_of_atom = self.__get_interaction_constants(constants_of_pairs, atom_i_type, atom_j_type)

        Vsss = integral_const_of_atom['V_sssigma']
        Vsps = integral_const_of_atom['V_spsigma']
        Vsds = integral_const_of_atom['V_sdsigma']
        Vstarsssigma = integral_const_of_atom['V_starsssigma']
        Vstarss = integral_const_of_atom['V_starssigma']
        Vstarps = integral_const_of_atom['V_starpsigma']
        Vstards = integral_const_of_atom['V_stardsigma']
        Vpps = integral_const_of_atom['V_ppsigma']
        Vppp = integral_const_of_atom['V_pppi']
        Vpds = integral_const_of_atom['V_pdsigma']
        Vpdp = integral_const_of_atom['V_pdpi']
        Vdds = integral_const_of_atom['V_ddsigma']
        Vddp = integral_const_of_atom['V_ddpi']
        Vddd = integral_const_of_atom['V_ddd']

        # ORDER OF MATRIX (little cheat sheet)
        # [0, 1,  2,  3,  4,   5,    6,    7,     8,  9]
        # [s, px, py, pz, dxy, dyz, dxz, dx2dy2, dz2, s*]

        nx, ny, nz = n['nx'], n['ny'], n['nz']

        H_SK[0][0] = Vsss
        H_SK[0][1] = nx * Vsps
        H_SK[0][2] = ny * Vsps
        H_SK[0][3] = nz * Vsps
        H_SK[0][4] = np.sqrt(3) * nx * ny * Vsds
        H_SK[0][5] = np.sqrt(3) * ny * nz * Vsds
        H_SK[0][6] = np.sqrt(3) * nx * nz * Vsds
        H_SK[0][7] = np.sqrt(3) / 2. * (nx * nx - ny * ny) * Vsds
        H_SK[0][8] = - 1. / 2. * (nx * nx + ny * ny - 2 * nz * nz) * Vsds
        
        H_SK[1][1] = nx * nx * Vpps + (1 - nx * nx) * Vppp
        H_SK[1][2] = - nx * ny * (Vppp - Vpps)
        H_SK[1][3] = - nx * nz * (Vppp - Vpps)
        H_SK[1][4] = np.sqrt(3) * nx * nx * ny * Vpds + (1 - 2 * nx * nx) * ny * Vpdp
        H_SK[1][5] = nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
        H_SK[1][6] = np.sqrt(3) * nx * nx * nz * Vpds + (1 - 2 * nx * nx) * nz * Vpdp
        H_SK[1][7] = np.sqrt(3) * nx * (nx * nx - ny * ny) * Vpds + nx * (1 - nx * nx + ny * ny) * Vpdp
        H_SK[1][8] = - np.sqrt(3) * nx * nz * nz * Vpdp - 0.5 * nx * (nx * nx + ny * ny - 2 * nz * nz) * Vpds

        H_SK[2][2] = ny * ny * Vpps + (1 - ny * ny) * Vppp
        H_SK[2][3] = - ny * nz * (Vppp - Vpps)
        H_SK[2][4] = np.sqrt(3) * ny * ny * nx * Vpds + (1 - 2 * ny * ny) * nx * Vpdp
        H_SK[2][5] = np.sqrt(3) * ny * ny * nz * Vpds + (1 - 2 * ny * ny) * nz * Vpdp
        H_SK[2][6] = nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
        H_SK[2][7] = np.sqrt(3) / 2. * ny * (nx * nx - ny * ny) * Vpds - ny * (1 - ny * ny + nx * nx) * Vpdp
        H_SK[2][8] = - np.sqrt(3) * ny * nz * nz * Vpdp - 0.5 * ny * (nx * nx + ny * ny - 2 * nz * nz) * Vpds

        H_SK[3][3] = nz * nz * Vpps + (1 - nz * nz) * Vppp
        H_SK[3][4] = nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
        H_SK[3][5] = np.sqrt(3) * nz * nz * ny * Vpds + (1 - 2 * nz * nz) * ny * Vpdp
        H_SK[3][6] = np.sqrt(3) * nz * nz * nx * Vpds + (1 - 2 * nz * nz) * nx * Vpdp
        H_SK[3][7] = np.sqrt(3) / 2. * nz * (nx * nx - ny * ny) * Vpds - nz * (nx * nx - ny * ny) * Vpdp
        H_SK[3][8] = np.sqrt(3) * nz * (nx * nx + ny * ny) * Vpdp - 0.5 * nz * (nx * nx + ny * ny - 2 * nz * nz) * Vpds

        H_SK[4][4] = nx**2 * ny**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + ny**2) * Vddp + nz**2 * Vddd
        H_SK[5][5] = nz**2 * ny**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nz**2 + ny**2) * Vddp + nx**2 * Vddd
        H_SK[6][6] = nx**2 * nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + nz**2) * Vddp + ny**2 * Vddd
        H_SK[4][5] = ny**2 * nx * nz * (3 * Vdds - 4 * Vddp + Vddd) + nx * nz * (Vddp - Vddd)
        H_SK[4][6] = nx**2 * ny * nz * (3 * Vdds - 4 * Vddp + Vddd) + ny * nz * (Vddp - Vddd)
        H_SK[5][6] = nz**2 * nx * ny * (3 * Vdds - 4 * Vddp + Vddd) + nx * ny * (Vddp - Vddd)
        H_SK[6][7] = 0.5 * nx * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) + 2 * (Vddp - Vddd))
        H_SK[5][7] = 0.5 * ny * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) - 2 * (Vddp - Vddd))
        H_SK[4][7] = 0.5 * nx * ny * (nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd)
        H_SK[4][8] = 0.5 * np.sqrt(3) * nx * ny * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)
        H_SK[5][8] = - 0.5 * np.sqrt(3) * ny * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))
        H_SK[6][8] = - 0.5 * np.sqrt(3) * nx * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))

        H_SK[7][7] = 0.25 * (nx**2 - ny**2)**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + ny**2) * Vddp + nz**2 * Vdds
        H_SK[8][8] = 0.75 * (nx**2 + ny**2)**2 * Vddd + 3 * (nx**2 + ny**2) * nz**2 * Vddp + 0.25 * (nx**2 + ny**2 - 2*nz**2)**2 * Vdds
        H_SK[7][8] = 0.25 * (nx**2 - ny**2) * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)

        H_SK[0][9] = Vstarsssigma
        H_SK[9][1] = nx * Vstarps
        H_SK[9][2] = ny * Vstarps
        H_SK[9][3] = nz * Vstarps
        H_SK[9][4] = np.sqrt(3) * nx * ny * Vstards
        H_SK[9][5] = np.sqrt(3) * ny * nz * Vstards
        H_SK[9][6] = np.sqrt(3) * nx * nz * Vstards
        H_SK[6][9] = np.sqrt(3) * nx * nz * Vstards
        H_SK[9][7] = np.sqrt(3) / 2. * (nx * nx - ny * ny) * Vstards
        H_SK[9][8] = - 1. / 2. * (nx * nx + ny * ny - 2 * nz * nz) * Vstards
        H_SK[9][9] = Vstarss
        
        H_SK[1][0] = - H_SK[0][1]
        H_SK[2][0] = - H_SK[0][2]
        H_SK[3][0] = - H_SK[0][3]
        H_SK[4][0] = H_SK[0][4]
        H_SK[5][0] = H_SK[0][5]
        H_SK[6][0] = H_SK[0][6]
        H_SK[7][0] = H_SK[0][7]
        H_SK[8][0] = H_SK[0][8]
        H_SK[2][1] = H_SK[1][2]
        H_SK[3][1] = H_SK[1][3]
        H_SK[4][1] = - H_SK[1][4]
        H_SK[5][1] = H_SK[1][5]
        H_SK[6][1] = - H_SK[1][6]
        H_SK[7][1] = - H_SK[1][7]
        H_SK[8][1] = - H_SK[1][8]
        H_SK[3][2] = H_SK[2][3]
        H_SK[4][2] = - H_SK[2][4]
        H_SK[5][2] = - H_SK[2][5]
        H_SK[6][2] = - H_SK[2][6]
        H_SK[7][2] = - H_SK[2][7]
        H_SK[8][2] = - H_SK[2][8]
        H_SK[4][3] = - H_SK[3][4]
        H_SK[5][3] = - H_SK[3][5]
        H_SK[6][3] = - H_SK[3][6]
        H_SK[7][3] = - H_SK[3][7]
        H_SK[8][3] = - H_SK[3][8]
        H_SK[5][4] = H_SK[4][5]        
        H_SK[6][4] = H_SK[4][6] 
        H_SK[6][5] = H_SK[5][6]   
        H_SK[7][6] = H_SK[6][7]
        H_SK[7][5] = H_SK[5][7]  
        H_SK[7][4] = H_SK[4][7]
        H_SK[8][4] = H_SK[4][8]
        H_SK[8][5] = H_SK[5][8]
        H_SK[8][6] = H_SK[6][8]
        H_SK[8][7] = H_SK[7][8]
        H_SK[1][9] = - H_SK[9][1]
        H_SK[2][9] = - H_SK[9][2]
        H_SK[3][9] = - H_SK[9][3]
        H_SK[4][9] = H_SK[9][4]
        H_SK[5][9] = H_SK[9][5]
        H_SK[6][9] = H_SK[9][6]
        H_SK[7][9] = H_SK[9][7]
        H_SK[8][9] = H_SK[9][8]
        H_SK[9][0] = - H_SK[0][9]

        return np.matrix(H_SK)

    def __spin_orbit_hamiltonian(self, atom_store: dict, atom_type: str, lp: float, ld: float, sgn: (float, None),
                                 sigma: str) -> np.matrix:

        """
        Method calculates diagonal matrix of band energies for spin-orbit interactions
        Args:
            atom_store: dict with physical constants (energies of bands) of atoms in lattice
            atom_type: type of element (for example 'C')
            lp: p-band interaction constant for spin-orbit interactions
            ld: d-band interaction constant for spin-orbit interactions
            sgn: spin sign
            sigma: spin
        Returns: Slater-Koster matrix for spin-orbit interactions
        """

        H_SK = np.zeros((10, 10), complex)

        state_energies_of_atom = self.get_atom_type(atom_store, atom_type)

        if state_energies_of_atom is None:

            return np.matrix(H_SK)

        else:

            if sigma == ' up up':

                Es = state_energies_of_atom['Es up up']
                Epx = state_energies_of_atom['Epx up up']
                Epy = state_energies_of_atom['Epy up up']
                Epz = state_energies_of_atom['Epz up up']
                Edz2 = state_energies_of_atom['Edz2 up up']
                Edxz = state_energies_of_atom['Edxz up up']
                Edyz = state_energies_of_atom['Edyz up up']
                Edxy = state_energies_of_atom['Edxy up up']
                Edx2y2 = state_energies_of_atom['Edx2y2 up up']
                Estar = state_energies_of_atom['Estar up up']

                H_SK[0][0] = Es
                H_SK[1][1] = Epx * sgn * lp
                H_SK[1][2] = - 1j * sgn * lp * 0.5
                H_SK[1][3] = 0

                H_SK[2][1] = 1j * 0.5 * sgn * lp
                H_SK[2][2] = Epy * sgn * lp
                H_SK[2][3] = 0

                H_SK[3][1] = 0
                H_SK[3][2] = 0
                H_SK[3][3] = Epz * sgn * lp

                H_SK[4][4] = Edxy
                H_SK[4][5] = 0
                H_SK[4][6] = 0
                H_SK[4][7] = 0
                H_SK[4][8] = 0

                H_SK[5][4] = 0
                H_SK[5][5] = Edyz
                H_SK[5][6] = 0
                H_SK[5][7] = 1j * 0.5 * sgn * ld
                H_SK[5][8] = 0

                H_SK[6][4] = 0
                H_SK[6][5] = 0
                H_SK[6][6] = Edz2
                H_SK[6][7] = 0
                H_SK[6][8] = 0

                H_SK[7][4] = 0
                H_SK[7][5] = -1j * 0.5 * sgn * ld
                H_SK[7][6] = 0
                H_SK[7][7] = Edxz
                H_SK[7][8] = 0

                H_SK[8][4] = -1j * sgn * ld
                H_SK[8][5] = 0
                H_SK[8][6] = 0
                H_SK[8][7] = 0
                H_SK[8][8] = Edx2y2

                H_SK[9][9] = Estar

            elif sigma == ' down down':

                Es = state_energies_of_atom['Es down down']
                Epx = state_energies_of_atom['Epx down down']
                Epy = state_energies_of_atom['Epy down down']
                Epz = state_energies_of_atom['Epz down down']
                Edz2 = state_energies_of_atom['Edz2 down down']
                Edxz = state_energies_of_atom['Edxz down down']
                Edyz = state_energies_of_atom['Edyz down down']
                Edxy = state_energies_of_atom['Edxy down down']
                Edx2y2 = state_energies_of_atom['Edx2y2 down down']
                Estar = state_energies_of_atom['Estar down down']

                H_SK[0][0] = Es
                H_SK[1][1] = Epx * sgn * lp
                H_SK[1][2] = - 1j * sgn * lp * 0.5
                H_SK[1][3] = 0

                H_SK[2][1] = 1j * 0.5 * sgn * lp
                H_SK[2][2] = Epy * sgn * lp
                H_SK[2][3] = 0

                H_SK[3][1] = 0
                H_SK[3][2] = 0
                H_SK[3][3] = Epz * sgn * lp

                H_SK[4][4] = Edxy
                H_SK[4][5] = 0
                H_SK[4][6] = 0
                H_SK[4][7] = 0
                H_SK[4][8] = 0

                H_SK[5][4] = 0
                H_SK[5][5] = Edyz
                H_SK[5][6] = 0
                H_SK[5][7] = 1j * 0.5 * sgn * ld
                H_SK[5][8] = 0

                H_SK[6][4] = 0
                H_SK[6][5] = 0
                H_SK[6][6] = Edz2
                H_SK[6][7] = 0
                H_SK[6][8] = 0

                H_SK[7][4] = 0
                H_SK[7][5] = -1j * 0.5 * sgn * ld
                H_SK[7][6] = 0
                H_SK[7][7] = Edxz
                H_SK[7][8] = 0

                H_SK[8][4] = -1j * sgn * ld
                H_SK[8][5] = 0
                H_SK[8][6] = 0
                H_SK[8][7] = 0
                H_SK[8][8] = Edx2y2

                H_SK[9][9] = Estar

            else:

                H_SK[0][0] = 0
                H_SK[1][1] = 0
                H_SK[1][2] = 0
                H_SK[1][3] = 1/2 * lp

                H_SK[2][1] = 0
                H_SK[2][2] = 0
                H_SK[2][3] = -1j * 0.5 * lp

                H_SK[3][1] = 0
                H_SK[3][1] = -1j * 0.5 * lp
                H_SK[3][2] = 1j * 0.5 * lp

                H_SK[4][4] = 0
                H_SK[4][5] = 1/2 * ld
                H_SK[4][6] = 0
                H_SK[4][7] = -1j * 0.5 * ld
                H_SK[4][8] = 0

                H_SK[5][4] = -1/2 * ld
                H_SK[5][5] = 0
                H_SK[5][6] = -1j * np.sqrt(3) * ld / 2
                H_SK[5][7] = 0
                H_SK[5][8] = -1j * 0.5 * ld

                H_SK[6][4] = 1j * 0.5 * ld
                H_SK[6][5] = 0
                H_SK[6][6] = np.sqrt(3) / 2 * ld
                H_SK[6][7] = 0
                H_SK[6][8] = -1j * 0.5 * ld

                H_SK[7][4] = 1j * 0.5 * ld
                H_SK[7][5] = 0
                H_SK[7][6] = np.sqrt(3) / 2 * ld
                H_SK[7][7] = 0
                H_SK[7][8] = -1j * 0.5 * ld

                H_SK[8][4] = 0
                H_SK[8][5] = 1j * 0.5 * ld
                H_SK[8][6] = 0
                H_SK[8][7] = 1/2 * ld
                H_SK[8][8] = 0

                H_SK[9][9] = 0

            new_order = [0, 1, 2, 3, 4, 5, 7, 8, 6, 9]
            H_SK = H_SK[:, new_order]
            
            H_SK[[6, 7]] = H_SK[[7, 8]]
            H_SK[[6, 8]] = H_SK[[8, 6]]
            
            return np.matrix(H_SK)

    def __calculate_energy_matrix(self, atom_store: dict, atom_type: str)->np.array:

        """
        Method calculates diagonal matrix of energies for band structure
        Args:
            atom_store: dict with physical constants (energies of bands) of atoms in lattice
            atom_type: type of element (for example 'C')
        Returns: diagonal matrix of band energies
        """

        atom_parameters = self.get_atom_type(atom_store, atom_type)
        H_Rii = np.zeros((10, 10))

        if atom_parameters is None:
            np.fill_diagonal(H_Rii, 0)

        else:
            atom_parameters = self.get_atom_type(atom_store, atom_type)

            Es = atom_parameters['Es']
            Epx = atom_parameters['Epx']
            Epy = atom_parameters['Epy']
            Epz = atom_parameters['Epz']
            Edz2 = atom_parameters['Edz2']
            Edxz = atom_parameters['Edxz']
            Edyz = atom_parameters['Edyz']
            Edxy = atom_parameters['Edxy']
            Edx2y2 = atom_parameters['Edx2y2']
            Estar = atom_parameters['Estar']

            row, col = np.diag_indices_from(H_Rii)
            H_Rii[row, col] = np.array([[Es, Epx, Epy, Epz, Edxy, Edyz, Edxz, Edx2y2, Edz2, Estar]])

        return H_Rii

    def calculate_spin_mixing_sk(self, calculation_type: str, ri: np.array,
                                 rj: np.array, constants_of_pairs: dict, atom_i: str, atom_j: str,
                                 flat: bool)->np.matrix:

        """
        Method returns interaction matrix of required shape. If calculation_type is 'non spin' then Slater Koster matrix
            is 10x10, else 20x20 (because of spin-orbit interactions)
        Args:
            calculation_type: If calculation_type is 'non spin' then Slater Koster matrix is 10x10, else 20x20.
            ri: position of i-th atom in lattice
            rj: position of i-th atom in lattice
            constants_of_pairs: dictionary with interaction constants
            atom_i: type of element of i-th atom in lattice
            atom_j: type of element of j-th atom in lattice
            flat: bool
        Returns: Slater Koster matrix of required dimension
        """

        if calculation_type == 'non spin':

            the_chosen_one = np.matrix(self.__calculate_slayterkoster_matrix(ri, rj, constants_of_pairs, atom_i, atom_j))
        else:

            equal_spin_matrix = self.__calculate_slayterkoster_matrix(ri, rj, constants_of_pairs, atom_i, atom_j)

            the_chosen_one = np.bmat([[equal_spin_matrix if i != j else np.zeros((10, 10))
                                       for i in range(2)] for j in range(2)])

        return the_chosen_one

    def calculate_spin_mixing_diagonal(self, calculation_type: str, atom_store: dict, atom_type: str, lp: float,
                                       ld: float) -> np.matrix:

        """
        Method returns matrix of diagonal energies of required shape. If calculation_type is 'non spin' then
            diagonal matrix is 10x10, else 20x20 (because of spin-orbit interactions).
        Args:
            calculation_type: If calculation_type is 'non spin' then diagonal matrix with band energies
            is 10x10, else 20x20 (because of spin-orbit interactions)
            atom_store: dict with physical constants (energies of bands) of atoms in lattice
            atom_type: type of element (for example 'C')
            lp: p-band interaction constant for spin-orbit interactions
            ld: d-band interaction constant for spin-orbit interactions
        Returns:
        """

        if calculation_type == 'non spin':

            diagonal_one = self.__calculate_energy_matrix(atom_store, atom_type)

        else:

            diagonal_one = np.bmat([[self.__spin_orbit_hamiltonian(atom_store, atom_type, lp, ld, 1, ' up up'),
                                     self.__spin_orbit_hamiltonian(atom_store, atom_type, lp, ld, None, ' up down')],
                                    [self.__spin_orbit_hamiltonian(atom_store, atom_type, lp, ld, -1, ' down down'),
                                     self.__spin_orbit_hamiltonian(atom_store, atom_type, lp, ld, None, ' up down').H]])

        return np.matrix(diagonal_one)

