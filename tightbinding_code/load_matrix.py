import numpy as np 
def __calculate_slayterkoster_matrix(self, ri: int, rj: int,
                                     constants_of_pairs: dict, atom_i_type: str, atom_j_type: str, flat: bool) -> np.matrix:

    """
    Method calculates Slater-Koster matrix (10x10)
    Args:
        ri: position of i-th atom in lattice
        rj: position of j-th atom in lattice
        constants_of_pairs: dictionary with interaction constants
        atom_i_type: type of element i-th atom
        atom_j_type: type of element i-th atom
        flat: bool type;
    Returns: Slater-Koster matrix
    """

    H_SK = np.zeros((10, 10))
    n = self.__get_directional_cosines(ri, rj)
    integral_const_of_atom = self.__get_interaction_constants(constants_of_pairs, atom_i_type, atom_j_type)

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



def __calculate_slayterkoster_matrix2(constants_of_pairs):

    """
    Method calculates Slater-Koster matrix (10x10)
    Args:
        ri: position of i-th atom in lattice
        rj: position of j-th atom in lattice
        constants_of_pairs: dictionary with interaction constants
        atom_i_type: type of element i-th atom
        atom_j_type: type of element i-th atom
        flat: bool type;
    Returns: Slater-Koster matrix
    """

    H_SK = np.zeros((10, 10))
    #n = self.__get_directional_cosines(ri, rj)
    integral_const_of_atom = constants_of_pairs

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

    nx, ny, nz = 1,1,1
    
    H_SK[0][0] = Vsss
    H_SK[0][1] = nx * Vsps
    H_SK[0][2] = ny * Vsps
    H_SK[0][3] = nz * Vsps
    H_SK[1][0] = - nx * Vsps
    H_SK[2][0] = - ny * Vsps
    H_SK[3][0] = - nz * Vsps
    H_SK[0][4] = np.sqrt(3) * nx * ny * Vsds
    H_SK[0][5] = np.sqrt(3) * ny * nz * Vsds
    H_SK[0][6] = np.sqrt(3) * nx * nz * Vsds
    H_SK[4][0] = np.sqrt(3) * nx * ny * Vsds
    H_SK[5][0] = np.sqrt(3) * ny * nz * Vsds
    H_SK[6][0] = np.sqrt(3) * nx * nz * Vsds
    H_SK[0][7] = np.sqrt(3) / 2. * (nx * nx - ny * ny) * Vsds
    H_SK[7][0] = np.sqrt(3) / 2. * (nx * nx - ny * ny) * Vsds
    H_SK[0][8] = - 1. / 2. * (nx * nx + ny * ny - 2 * nz * nz) * Vsds
    H_SK[8][0] = - 1. / 2. * (nx * nx + ny * ny - 2 * nz * nz) * Vsds

    H_SK[1][1] = nx * nx * Vpps + (1 - nx * nx) * Vppp
    H_SK[1][2] = - nx * ny * (Vppp - Vpps)
    H_SK[2][1] = - nx * ny * (Vppp - Vpps)
    H_SK[1][3] = - nx * nz * (Vppp - Vpps)
    H_SK[3][1] = - nx * nz * (Vppp - Vpps)
    H_SK[1][4] = np.sqrt(3) * nx * nx * ny * Vpds + (1 - 2 * nx * nx) * ny * Vpdp
    H_SK[1][5] = nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
    H_SK[1][6] = np.sqrt(3) * nx * nx * nz * Vpds + (1 - 2 * nx * nx) * nz * Vpdp
    H_SK[4][1] = - np.sqrt(3) * nx * nx * ny * Vpds - (1 - 2 * nx * nx) * ny * Vpdp
    H_SK[5][1] = nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
    H_SK[6][1] = - np.sqrt(3) * nx * nx * nz * Vpds - (1 - 2 * nx * nx) * nz * Vpdp
    H_SK[1][7] = np.sqrt(3) * nx * (nx * nx - ny * ny) * Vpds + nx * (1 - nx * nx + ny * ny) * Vpdp
    H_SK[7][1] = - np.sqrt(3) * nx * (nx * nx - ny * ny) * Vpds - nx * (1 - nx * nx + ny * ny) * Vpdp
    H_SK[1][8] = - np.sqrt(3) * nx * nz * nz * Vpdp - 0.5 * nx * (nx * nx + ny * ny - 2 * nz * nz) * Vpds
    H_SK[8][1] = np.sqrt(3) * nx * nz * nz * Vpdp + 0.5 * nx * (nx * nx + ny * ny - 2 * nz * nz) * Vpds

    H_SK[2][2] = ny * ny * Vpps + (1 - ny * ny) * Vppp
    H_SK[2][3] = - ny * nz * (Vppp - Vpps)
    H_SK[3][2] = - ny * nz * (Vppp - Vpps)
    H_SK[2][4] = np.sqrt(3) * ny * ny * nx * Vpds + (1 - 2 * ny * ny) * nx * Vpdp
    H_SK[4][2] = - np.sqrt(3) * ny * ny * nx * Vpds - (1 - 2 * ny * ny) * nx * Vpdp
    H_SK[2][5] = np.sqrt(3) * ny * ny * nz * Vpds + (1 - 2 * ny * ny) * nz * Vpdp
    H_SK[5][2] = - np.sqrt(3) * ny * ny * nz * Vpds - (1 - 2 * ny * ny) * nz * Vpdp
    H_SK[2][6] = nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
    H_SK[6][2] = - nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
    H_SK[2][7] = np.sqrt(3) / 2. * ny * (nx * nx - ny * ny) * Vpds - ny * (1 - ny * ny + nx * nx) * Vpdp
    H_SK[7][2] = - np.sqrt(3) / 2. * ny * (nx * nx - ny * ny) * Vpds + ny * (1 - ny * ny + nx * nx) * Vpdp
    H_SK[2][8] = - np.sqrt(3) * ny * nz * nz * Vpdp - 0.5 * ny * (nx * nx + ny * ny - 2 * nz * nz) * Vpds
    H_SK[8][2] = np.sqrt(3) * ny * nz * nz * Vpdp + 0.5 * ny * (nx * nx + ny * ny - 2 * nz * nz) * Vpds

    H_SK[3][3] = nz * nz * Vpps + (1 - nz * nz) * Vppp
    H_SK[3][4] = nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
    H_SK[4][3] = - nx * ny * nz * (np.sqrt(3) * Vpds - 2 * Vpdp)
    H_SK[3][5] = np.sqrt(3) * nz * nz * ny * Vpds + (1 - 2 * nz * nz) * ny * Vpdp
    H_SK[5][3] = - np.sqrt(3) * nz * nz * ny * Vpds - (1 - 2 * nz * nz) * ny * Vpdp
    H_SK[3][6] = np.sqrt(3) * nz * nz * nx * Vpds + (1 - 2 * nz * nz) * nx * Vpdp
    H_SK[6][3] = - np.sqrt(3) * nz * nz * nx * Vpds - (1 - 2 * nz * nz) * nx * Vpdp
    H_SK[3][7] = np.sqrt(3) / 2. * nz * (nx * nx - ny * ny) * Vpds - nz * (nx * nx - ny * ny) * Vpdp
    H_SK[7][3] = - np.sqrt(3) / 2. * nz * (nx * nx - ny * ny) * Vpds + nz * (nx * nx - ny * ny) * Vpdp
    H_SK[3][8] = np.sqrt(3) * nz * (nx * nx + ny * ny) * Vpdp - 0.5 * nz * (nx * nx + ny * ny - 2 * nz * nz) * Vpds
    H_SK[8][3] = - np.sqrt(3) * nz * (nx * nx + ny * ny) * Vpdp + 0.5 * nz * (nx * nx + ny * ny - 2 * nz * nz) * Vpds

    H_SK[4][4] = nx**2 * ny**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + ny**2) * Vddp + nz**2 * Vddd
    H_SK[5][5] = nz**2 * ny**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nz**2 + ny**2) * Vddp + nx**2 * Vddd
    H_SK[6][6] = nx**2 * nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + nz**2) * Vddp + ny**2 * Vddd
    H_SK[4][5] = ny**2 * nx * nz * (3 * Vdds - 4 * Vddp + Vddd) + nx * nz * (Vddp - Vddd)
    H_SK[5][4] = ny**2 * nx * nz * (3 * Vdds - 4 * Vddp + Vddd) + nx * nz * (Vddp - Vddd)
    H_SK[4][6] = nx**2 * ny * nz * (3 * Vdds - 4 * Vddp + Vddd) + ny * nz * (Vddp - Vddd)
    H_SK[6][4] = nx**2 * ny * nz * (3 * Vdds - 4 * Vddp + Vddd) + ny * nz * (Vddp - Vddd)
    H_SK[5][6] = nz**2 * nx * ny * (3 * Vdds - 4 * Vddp + Vddd) + nx * ny * (Vddp - Vddd)
    H_SK[6][5] = nz**2 * nx * ny * (3 * Vdds - 4 * Vddp + Vddd) + nx * ny * (Vddp - Vddd)
    H_SK[6][7] = 0.5 * nx * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) + 2 * (Vddp - Vddd))
    H_SK[7][6] = 0.5 * nx * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) + 2 * (Vddp - Vddd))
    H_SK[5][7] = 0.5 * ny * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) - 2 * (Vddp - Vddd))
    H_SK[7][5] = 0.5 * ny * nz *((nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd) - 2 * (Vddp - Vddd))
    H_SK[4][7] = 0.5 * nx * ny * (nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd)
    H_SK[7][4] = 0.5 * nx * ny * (nx**2 - ny**2) * (3 * Vdds - 4 * Vddp + Vddd)
    H_SK[4][8] = 0.5 * np.sqrt(3) * nx * ny * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)
    H_SK[5][8] = - 0.5 * np.sqrt(3) * ny * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))
    H_SK[6][8] = - 0.5 * np.sqrt(3) * nx * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))
    H_SK[8][4] = 0.5 * np.sqrt(3) * nx * ny * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)
    H_SK[8][5] = - 0.5 * np.sqrt(3) * ny * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))
    H_SK[8][6] = - 0.5 * np.sqrt(3) * nx * nz * ((nx**2 + ny**2) * (Vdds - 2 * Vddp + Vdds) + 2 * nz**2 * (Vddp - Vdds))

    H_SK[7][7] = 0.25 * (nx**2 - ny**2)**2 * (3 * Vdds - 4 * Vddp + Vddd) + (nx**2 + ny**2) * Vddp + nz**2 * Vdds
    H_SK[8][8] = 0.75 * (nx**2 + ny**2)**2 * Vddd + 3 * (nx**2 + ny**2) * nz**2 * Vddp + 0.25 * (nx**2 + ny**2 - 2*nz**2)**2 * Vdds
    H_SK[7][8] = 0.25 * (nx**2 - ny**2) * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)
    H_SK[8][7] = 0.25 * (nx**2 - ny**2) * (nz**2 * (3 * Vdds - 4 * Vddp + Vddd) + Vddd - Vdds)

    H_SK[0][9] = Vstarss * Vsss
    H_SK[9][0] = -Vstarss * Vsss
    H_SK[9][1] = nx * Vstarps
    H_SK[9][2] = ny * Vstarps
    H_SK[9][3] = nz * Vstarps
    H_SK[1][9] = - nx * Vstarps
    H_SK[2][9] = - ny * Vstarps
    H_SK[3][9] = - nz * Vstarps
    H_SK[9][4] = np.sqrt(3) * nx * ny * Vstards
    H_SK[9][5] = np.sqrt(3) * ny * nz * Vstards
    H_SK[9][6] = np.sqrt(3) * nx * nz * Vstards
    H_SK[4][9] = np.sqrt(3) * nx * ny * Vstards
    H_SK[5][9] = np.sqrt(3) * ny * nz * Vstards
    H_SK[6][9] = np.sqrt(3) * nx * nz * Vstards
    H_SK[9][7] = np.sqrt(3) / 2. * (nx * nx - ny * ny) * Vstards
    H_SK[7][9] = np.sqrt(3) / 2. * (nx * nx - ny * ny) * Vstards
    H_SK[9][8] = - 1. / 2. * (nx * nx + ny * ny - 2 * nz * nz) * Vstards
    H_SK[8][9] = - 1. / 2. * (nx * nx + ny * ny - 2 * nz * nz) * Vstards
    H_SK[9][9] = Vstarss

    return np.matrix(H_SK)
    
    
def __calculate_slayterkoster_matrix1(constants_of_pairs):

    """
    Method calculates Slater-Koster matrix (10x10)
    Args:
        ri: position of i-th atom in lattice
        rj: position of j-th atom in lattice
        constants_of_pairs: dictionary with interaction constants
        atom_i_type: type of element i-th atom
        atom_j_type: type of element i-th atom
        flat: bool type;
    Returns: Slater-Koster matrix
    """

    H_SK = np.zeros((10, 10))
    #n = self.__get_directional_cosines(ri, rj)
    integral_const_of_atom = constants_of_pairs

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

    nx, ny, nz = 1,1,1
    
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

    H_SK[0][9] = Vstarss * Vsss
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
    
constans_of_pairs_example = {'V_sssigma': 1.,
                                       'V_spsigma': 2.,
                                       'V_sdsigma': 3.,
                                       'V_starssigma': 4.,
                                       'V_starpsigma': 5.,
                                       'V_stardsigma': 6.,
                                       'V_ppsigma': 7.,
                                       'V_pppi': -2.6,
                                       'V_pdsigma': 8.,
                                       'V_pdpi': 9.,
                                       'V_ddsigma': 10.,
                                       'V_ddpi': 11.,
                                       'V_ddd': 12.}
                                       
a = __calculate_slayterkoster_matrix1(constans_of_pairs_example)
b = __calculate_slayterkoster_matrix2(constans_of_pairs_example)
c = a == b
print(c)
