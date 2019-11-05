diagonal_energies = {'C': {'Es': -2.990,
                           'Epx': 3.710,
                           'Epy': 3.710,
                           'Epz': 3.710,
                           'Edz2': 400.,
                           'Edxz': 400.,
                           'Edyz': 400.,
                           'Edxy': 400.,
                           'Edx2y2': 400.,
                           'Estar': 400.}}   
                                           

interactive_constans = {('C', 'C'): {'V_sssigma': -5.000,
                                     'V_spsigma': 4.700,
                                     'V_sdsigma': 0.,
                                     'V_starssigma': 0.,
                                     'V_starpsigma': 0.,
                                     'V_stardsigma': 0.,
                                     'V_ppsigma': 5.500,
                                     'V_pppi': -1.550,
                                     'V_pdsigma': 0.,
                                     'V_pdpi': 0.,
                                     'V_ddsigma': 0.,
                                     'V_ddpi': 0.,
                                     'V_ddd': 0.}
                                     } 


interactive_constans_spin_spin_orbita = {'C': {'Es up up': 1,
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

configuration = {'parametrization':
{'lattice_type': 'rectangular_new',
 'gauss_sigma': 0.015,
 'start': -10,
 'stop': 10,
 'step': 0.001,
 'distance': float(format(1.536329, '.12g')),
 'lanczos_vectors':None,
 'vertical_num_of_steps': 20,
 'horizontal_num_of_steps': 70,
 'x_num_of_steps': 39,
 'sigma': 0.0001,
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 3000,
 'interactive_constans': interactive_constans,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': None,
 'number_of_defects': None
}

} 


# TODO policzyc uklady z parametryzacją publikacje 
# TODO poiczyć dla dewastacji na brzegach heksagonalnej struktury i prostokątnej
