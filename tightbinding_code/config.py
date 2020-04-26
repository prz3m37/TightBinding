diagonal_energies = {'C': {'Es': -50,
                           'Epx': -50,
                           'Epy': -50,
                           'Epz': 0,
                           'Edz2': 4000.,
                           'Edxz': 4000.,
                           'Edyz': 4000.,
                           'Edxy': 4000.,
                           'Edx2y2': 4000.,
                           'Estar': 4000.}}


interactive_constants = {('C', 'C'): {'V_sssigma': 0,
                                      'V_spsigma': 0,
                                      'V_sdsigma': 0.,
                                      'V_starsssigma': 0,
                                      'V_starssigma': 0.,
                                      'V_starpsigma': 0.,
                                      'V_stardsigma': 0.,
                                      'V_ppsigma': 0,
                                      'V_pppi': -2.6,
                                      'V_pdsigma': 0.,
                                      'V_pdpi': 0.,
                                      'V_ddsigma': 0.,
                                      'V_ddpi': 0.,
                                      'V_ddd': 0.}}

settings = {
    'data_source': None,
    'data_identification': None,
    'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
    'input data': "",
    'data base': None,
    'host': None,
    'port': None,
    'password': None,
    'db table': None,
    'select': None}


configuration = {
    'parametrization': {
        'ld': None,
        'lp': None,
        'power': None,
        'epsilon': None,
        'magnitude': 'LM',
        "number_of_cpus": 1,
        'spin_included': False,
        'distance': 1.42,
        'lanczos_vectors': None,
        'number_of_eigenvalues': 1,
        'on_site_energies': diagonal_energies,
        'interactive_constants': interactive_constants,
        'fermi_level': 0}
}
