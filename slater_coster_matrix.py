import numpy as np


class SlaterKoster(object):

    """
    SlaterKoster class calculates Slater Koster matrix which describes interactions between two atoms in lattice.
        User can define whether interactions wants to take Spin Orbit interaction into consideration.

    """

    @staticmethod
    def get_atom_type(atom_store: dict, atom_type: str):

        """

        Args:
            atom_store:
            atom_type:

        Returns:

        """

        if atom_store is None:
            atom = None
        else:
            atom = atom_store[atom_type]

        return atom

    @staticmethod
    def __get_interaction_constans(constants: dict, atom_i_type: str, atom_j_type: str):

        """

        Args:
            constants:
            atom_i_type:
            atom_j_type:

        Returns:

        """

        if constants is None:

            const = {'V_sssigma': 0,
                     'V_spsigma': 0,
                     'V_sdsigma': 0,
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
    def __get_directional_cosines(ri: np.array, rj: np.array):

        """

        Args:
            ri:
            rj:

        Returns:

        """

        r_diff = ri - rj
        r_diff_mod = np.linalg.norm(ri - rj)

        n = {'nx': r_diff[0] / r_diff_mod, 'ny': r_diff[1] / r_diff_mod, 'nz': r_diff[2] / r_diff_mod}
        return n

    def __calculate_slayterkoster_matrix(self, ri, rj, constants_of_pairs, atom_i_type, atom_j_type, flat):

        """

        Args:
            ri:
            rj:
            constants_of_pairs:
            atom_i_type:
            atom_j_type:
            flat:

        Returns:

        """

        H_SK = np.zeros((10, 10))
        n = self.__get_directional_cosines(ri, rj)
        integral_const_of_atom = self.__get_interaction_constans(constants_of_pairs, atom_i_type, atom_j_type)

        Vsss = integral_const_of_atom['V_sssigma']
        Vsps = integral_const_of_atom['V_spsigma']
        Vsds = integral_const_of_atom['V_sdsigma']
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

        # [0, 1,  2,  3,  4,   5,    6,    7,     8,  9]
        # [s, px, py, pz, dxy, dyz, dxz, dx2dy2, dz2, s*] old
        # [s, px, py, pz, dxy, dyz, dz2, dxz, dx2dy2, s*] new

        nx, ny, nz = n['nx'], n['ny'], n['nz']

        H_SK[0][0] = Vsss
        H_SK[0][1] = nx * Vsps
        H_SK[0][2] = ny * Vsps
        H_SK[0][3] = nz * Vsps
        H_SK[0][4] = np.sqrt(3) * nx * ny * Vsds
        H_SK[0][5] = np.sqrt(3) * ny * nz * Vsds
        #     H_SK[0][6] = np.sqrt(3) * nx * nz * Vsds
        #     H_SK[0][7] = np.sqrt(3) / 2. * (nx**2 - ny**2) * Vsds
        #     H_SK[0][8] = - 1. / 2. * (nx**2 + ny**2 - 2 * nz**2) * Vsds

        H_SK[0][8] = np.sqrt(3) * nx * nz * Vsds
        H_SK[0][6] = np.sqrt(3) / 2. * (nx**2 - ny**2) * Vsds
        H_SK[0][7] = - 1. / 2. * (nx**2 + ny**2 - 2 * nz**2) * Vsds

        H_SK[0][9] = Vstarss * Vsss

        H_SK[1][9] = nx * Vstarps
        H_SK[2][9] = ny * Vstarps
        H_SK[3][9] = nz * Vstarps
        H_SK[4][9] = np.sqrt(3) * nx * ny * Vstards
        H_SK[5][9] = np.sqrt(3) * ny * nz * Vstards
        #     H_SK[6][9] = np.sqrt(3) * nx * nz * Vstards
        #     H_SK[7][9] = np.sqrt(3) / 2 * (nx**2 - ny**2) * Vstards
        #     H_SK[8][9] = - 1 / 2 * (nx**2 + ny**2 - 2 * nz**2) * Vstards

        H_SK[8][9] = np.sqrt(3) * nx * nz * Vstards
        H_SK[6][9] = np.sqrt(3) / 2 * (nx**2 - ny**2) * Vstards
        H_SK[7][9] = - 1 / 2 * (nx**2 + ny**2 - 2 * nz**2) * Vstards

        H_SK[9][9] = Vstarss

        if flat is False:
            H_SK[1][1] = nx**2 * Vpps + (1 - nx**2) * Vppp
            H_SK[1][2] = - nx * ny * (Vppp - Vpps)
            H_SK[2][2] = ny**2 * Vpps + (1 - ny**2) * Vppp
        else:
            H_SK[1][1] = 0
            H_SK[1][2] = 0
            H_SK[2][2] = 0

        H_SK[1][3] = - nx * nz * (Vppp - Vpps)
        H_SK[1][4] = np.sqrt(3) * nx**2 * ny * Vpds + (1 - 2 * nx**2) * ny * Vpdp
        H_SK[1][5] = nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
        #     H_SK[1][6] = np.sqrt(3) * nx**2 * nz * Vpds + (1 - 2 * nx**2) * nz * Vpdp
        #     H_SK[1][7] = np.sqrt(3) * nx * (nx**2 - ny**2) * Vpds + nx * (1 - nx**2 + ny**2) * Vpdp
        #     H_SK[1][8] = - np.sqrt(3) * nx * nz**2 * Vpdp - 0.5 * nx * (nx**2 + ny**2 - 2 * nz**2) * Vpds

        H_SK[1][8] = np.sqrt(3) * nx**2 * nz * Vpds + (1 - 2 * nx**2) * nz * Vpdp
        H_SK[1][6] = np.sqrt(3) * nx * (nx**2 - ny**2) * Vpds + nx * (1 - nx**2 + ny**2) * Vpdp
        H_SK[1][7] = - np.sqrt(3) * nx * nz**2 * Vpdp - 0.5 * nx * (nx**2 + ny**2 - 2 * nz**2) * Vpds

        H_SK[2][3] = - ny * nz * (Vppp - Vpps)
        H_SK[2][4] = np.sqrt(3) * ny**2 * nx * Vpds + (1 - 2 * ny**2) * nx * Vpdp
        H_SK[2][5] = np.sqrt(3) * ny**2 * nz * Vpds + (1 - 2 * ny**2) * nz * Vpdp
        #     H_SK[2][6] = nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
        #     H_SK[2][7] = np.sqrt(3) / 2. * ny * (nx**2 - ny**2) * Vpds - ny * (1 - ny**2 + nx**2) * Vpdp
        #     H_SK[2][8] = - np.sqrt(3) * ny * nz**2 * Vpdp - 0.5 * ny * (nx**2 + ny**2 - 2 * nz**2) * Vpds

        H_SK[2][8] = nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
        H_SK[2][6] = np.sqrt(3) / 2. * ny * (nx**2 - ny**2) * Vpds - ny * (1 - ny**2 + nx**2) * Vpdp
        H_SK[2][7] = - np.sqrt(3) * ny * nz**2 * Vpdp - 0.5 * ny * (nx**2 + ny**2 - 2 * nz**2) * Vpds

        H_SK[3][3] = nz**2 * Vpps + (1 - nz**2) * Vppp
        H_SK[3][4] = nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
        H_SK[3][5] = np.sqrt(3) * nz**2 * ny * Vpds + (1 - 2 * nz**2) * ny * Vpdp
        #     H_SK[3][6] = np.sqrt(3) * nz**2 * nx * Vpds + (1 - 2 * nz**2) * nx * Vpdp
        #     H_SK[3][7] = np.sqrt(3) / 2. * nz * (nx**2 - ny**2) * Vpds - nz * (nx**2 - ny**2) * Vpdp
        #     H_SK[3][8] = np.sqrt(3) * nz * (nx**2 + ny**2) * Vpdp - 0.5 * nz * (nx**2 + ny**2 - 2 * nz**2) * Vpds

        H_SK[3][8] = np.sqrt(3) * nz**2 * nx * Vpds + (1 - 2 * nz**2) * nx * Vpdp
        H_SK[3][6] = np.sqrt(3) / 2. * nz * (nx**2 - ny**2) * Vpds - nz * (nx**2 - ny**2) * Vpdp
        H_SK[3][7] = np.sqrt(3) * nz * (nx**2 + ny**2) * Vpdp - 0.5 * nz * (nx**2 + ny**2 - 2 * nz**2) * Vpds

        H_SK[4][4] = nx**2 * ny**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + ny**2) * Vddp + nz**2 * Vddd
        H_SK[5][5] = nz**2 * ny**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nz**2 + ny**2) * Vddp + nx**2 * Vddd

        #     H_SK[6][6] = nx**2 * nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + nz**2) * Vddp + ny**2 * Vddd
        #     H_SK[4][5] = ny**2 * nx * nz * (3 * Vdds - 4 * Vddp + Vddd) + nx * nz * (Vddp - Vddd)
        #     H_SK[4][6] = nx**2 * ny * nz * (3 * Vdds - 4 * Vddp + Vddd) + ny * nz * (Vddp - Vddd)
        #     H_SK[5][6] = nz**2 * nx * ny * (3 * Vdds - 4 * Vddp + Vddd) + nx * ny * (Vddp - Vddd)
        #     H_SK[6][7] = 0.5 * nx * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) + 2 * (Vddp - Vddd))
        #     H_SK[5][7] = 0.5 * ny * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) - 2 * (Vddp - Vddd))
        #     H_SK[4][7] = 0.5 * nx * ny * (nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd)
        #     H_SK[4][8] = 0.5 * np.sqrt(3) * nx * ny * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)
        #     H_SK[5][8] = - 0.5 * np.sqrt(3) * ny * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))
        #     H_SK[6][8] = - 0.5 * np.sqrt(3) * nx * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))

        H_SK[8][8] = nx**2 * nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + nz**2) * Vddp + ny**2 * Vddd
        H_SK[4][5] = ny**2 * nx * nz * (3 * Vdds - 4 * Vddp + Vddd) + nx * nz * (Vddp - Vddd)
        H_SK[4][8] = nx**2 * ny * nz * (3 * Vdds - 4 * Vddp + Vddd) + ny * nz * (Vddp - Vddd)
        H_SK[5][8] = nz**2 * nx * ny * (3 * Vdds - 4 * Vddp + Vddd) + nx * ny * (Vddp - Vddd)
        H_SK[8][6] = 0.5 * nx * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) + 2 * (Vddp - Vddd))
        H_SK[5][6] = 0.5 * ny * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) - 2 * (Vddp - Vddd))
        H_SK[4][6] = 0.5 * nx * ny * (nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd)
        H_SK[4][7] = 0.5 * np.sqrt(3) * nx * ny * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)
        H_SK[5][7] = - 0.5 * np.sqrt(3) * ny * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))
        H_SK[8][7] = - 0.5 * np.sqrt(3) * nx * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))

        #     H_SK[7][7] = 0.25 * (nx**2 - ny**2)**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + ny**2) * Vddp + nz**2 * Vdds
        #     H_SK[8][8] = 0.75 * (nx**2 + ny**2)**2 * Vddd + 3 * (nx**2 + ny**2) * nz**2 * Vddp + 0.25 * (nx**2 + ny**2 - 2*nz**2)**2 * Vdds
        H_SK[6][6] = 0.25 * (nx**2 - ny**2)**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + ny**2) * Vddp + nz**2 * Vdds
        H_SK[7][7] = 0.75 * (nx**2 + ny**2)**2 * Vddd + 3 * (nx**2 + ny**2) * nz**2 * Vddp + 0.25 * (nx**2 + ny**2 - 2*nz**2)**2 * Vdds

        H_SK[1][0] = - H_SK[0][1]
        H_SK[2][0] = - H_SK[0][2]
        H_SK[3][0] = - H_SK[0][3]
        H_SK[4][1] = - H_SK[1][4]
        H_SK[4][2] = - H_SK[2][4]
        H_SK[4][3] = - H_SK[3][4]
        H_SK[5][1] = - H_SK[1][5]
        H_SK[5][2] = - H_SK[2][5]
        H_SK[5][3] = - H_SK[3][5]
        H_SK[6][1] = - H_SK[1][6]
        H_SK[6][2] = - H_SK[2][6]
        H_SK[6][3] = - H_SK[3][6]
        H_SK[7][1] = - H_SK[1][7]
        H_SK[7][2] = - H_SK[2][7]
        H_SK[7][3] = - H_SK[3][7]
        H_SK[8][1] = - H_SK[1][8]
        H_SK[8][2] = - H_SK[2][8]
        H_SK[8][3] = - H_SK[3][8]
        H_SK[1][9] = - H_SK[9][1]
        H_SK[2][9] = - H_SK[9][2]
        H_SK[3][9] = - H_SK[9][3]

        H_SK[2][1] = H_SK[1][2]
        H_SK[3][1] = H_SK[1][3]
        H_SK[3][2] = H_SK[2][3]
        H_SK[4][0] = H_SK[0][4]
        H_SK[5][0] = H_SK[0][5]
        H_SK[6][0] = H_SK[0][6]
        H_SK[7][0] = H_SK[0][7]
        H_SK[8][0] = H_SK[0][8]
        H_SK[5][4] = H_SK[4][5]
        H_SK[6][4] = H_SK[4][6]
        H_SK[7][4] = H_SK[4][7]
        H_SK[8][4] = H_SK[4][8]
        H_SK[6][5] = H_SK[5][6]
        H_SK[7][5] = H_SK[5][7]
        H_SK[8][5] = H_SK[5][8]
        H_SK[7][6] = H_SK[6][7]
        H_SK[8][6] = H_SK[6][8]
        H_SK[8][7] = H_SK[7][8]
        H_SK[0][9] = H_SK[9][0]
        H_SK[4][9] = H_SK[9][4]
        H_SK[5][9] = H_SK[9][5]
        H_SK[6][9] = H_SK[9][6]
        H_SK[7][9] = H_SK[9][7]
        H_SK[8][9] = H_SK[9][8]

        return np.matrix(H_SK)

    def __spin_orbit_hamiltonian(self, atom_store, atom_type, lp, ld, sgn, sigma):

        """

        Args:
            atom_store:
            atom_type:
            lp:
            ld:
            sgn:
            sigma:

        Returns:

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

            return np.matrix(H_SK)

    def __calculate_energy_matrix(self, atom_store, atom_type):

        """

        Args:
            atom_store:
            atom_type:

        Returns:

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
            H_Rii[row, col] = np.array([[Es, Epx, Epy, Epz, Edxy, Edyz, Edz2, Edxz, Edx2y2, Estar]])

        return H_Rii

    def calculate_spin_mixing_sk(self, calculation_type, ri, rj, constants_of_pairs, atom_i, atom_j, flat):

        """

        Args:
            calculation_type:
            ri:
            rj:
            constants_of_pairs:
            atom_i:
            atom_j:
            flat:

        Returns:

        """

        if calculation_type == 'non spin':

            the_chosen_one = np.matrix(self.__calculate_slayterkoster_matrix(ri, rj, constants_of_pairs, atom_i, atom_j,
                                                                             flat))
        else:

            equal_spin_matrix = self.__calculate_slayterkoster_matrix(ri, rj, constants_of_pairs, atom_i, atom_j, flat)

            the_chosen_one = np.bmat([[equal_spin_matrix if i != j else np.zeros((10, 10))
                                       for i in range(2)] for j in range(2)])

        return the_chosen_one

    def calculate_spin_mixing_diagonal(self, calculation_type, atom_store, atom_type, lp, ld):

        """

        Args:
            calculation_type:
            atom_store:
            atom_type:
            lp:
            ld:

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

