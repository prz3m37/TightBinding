import numpy as np

def V_alpha(h, r=1.42):
  n = 2.
  n_c = 6.5
  r_o = 1.536329
  r_c = 2.18
  s = np.exp((-n*(r/r_c)**(n_c) + n*(r/r_c)**(n_c))) * (r_o/r)**n 
  h_alpha = h * s
  return h_alpha

diagonal_energies1 = {'C': {'Es': -50,
                           'Epx': -50,
                           'Epy': -50,
                           'Epz': 0,
                           'Edz2': 4000.,
                           'Edxz': 4000.,
                           'Edyz': 4000.,
                           'Edxy': 4000.,
                           'Edx2y2': 4000.,
                           'Estar': 4000.}}   

                              
interactive_constans1 = {('C', 'C'): {'V_sssigma': 0,
                                     'V_spsigma': 0,
                                     'V_sdsigma': 0.,
                                     'V_starsssigma':0,
                                     'V_starssigma': 0.,
                                     'V_starpsigma': 0.,
                                     'V_stardsigma': 0.,
                                     'V_ppsigma': 0,
                                     'V_pppi': -2.6,
                                     'V_pdsigma': 0.,
                                     'V_pdpi': 0.,
                                     'V_ddsigma': 0.,
                                     'V_ddpi': 0.,
                                     'V_ddd': 0.}
                                     } 


diagonal_energies2 = {'C': {'Es': -2.990,
                           'Epx': 3.710,
                           'Epy': 3.710,
                           'Epz': 3.710,
                           'Edz2': 4000.,
                           'Edxz': 4000.,
                           'Edyz': 4000.,
                           'Edxy': 4000.,
                           'Edx2y2': 4000.,
                           'Estar': 4000.}}   
                                           

interactive_constans2 = {('C', 'C'): {'V_sssigma': -5.000,
                                     'V_spsigma': 4.700,
                                     'V_sdsigma': 0.,
                                     'V_starsssigma':0,
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


diagonal_energies3 = {'C': {'Es': -10.290,
                           'Epx': 0,
                           'Epy': 0,
                           'Epz': 0,
                           'Edz2': 4000.,
                           'Edxz': 4000.,
                           'Edyz': 4000.,
                           'Edxy': 4000.,
                           'Edx2y2': 4000.,
                           'Estar': 4000.}}   
                                           

interactive_constans3 = {('C', 'C'): {'V_sssigma': -8.42256,
                                     'V_spsigma': 8.08162,
                                     'V_sdsigma': 0.,
                                     'V_starsssigma':0,
                                     'V_starssigma': 0.,
                                     'V_starpsigma': 0.,
                                     'V_stardsigma': 0.,
                                     'V_ppsigma': 7.75792,
                                     'V_pppi': -3.67510,
                                     'V_pdsigma': 0.,
                                     'V_pdpi': 0.,
                                     'V_ddsigma': 0.,
                                     'V_ddpi': 0.,
                                     'V_ddd': 0.}
                                     } 


diagonal_energies4 = {'C': {'Es': -2.99,
                           'Epx': 3.71,
                           'Epy': 3.71,
                           'Epz': 3.71,
                           'Edz2': 4000.,
                           'Edxz': 4000.,
                           'Edyz': 4000.,
                           'Edxy': 4000.,
                           'Edx2y2': 4000.,
                           'Estar': 4000.}}   
                                           

interactive_constans4 = {('C', 'C'): {'V_sssigma': V_alpha(h=-5., r=1.42),
                                     'V_spsigma': V_alpha(h=4.7, r=1.42),
                                     'V_sdsigma': 0.,
                                     'V_starsssigma':0,
                                     'V_starssigma': 0.,
                                     'V_starpsigma': 0.,
                                     'V_stardsigma': 0.,
                                     'V_ppsigma': V_alpha(h=5.5, r=1.42),
                                     'V_pppi': V_alpha(h=-1.55, r=1.42),
                                     'V_pdsigma': 0.,
                                     'V_pdpi': 0.,
                                     'V_ddsigma': 0.,
                                     'V_ddpi': 0.,
                                     'V_ddd': 0.}
                                     } 




diagonal_energies5 = {'C': {'Es': -1.0458,
                           'Epx': 7.0850,
                           'Epy': 7.0850,
                           'Epz': 7.0850,
                           'Edz2': 27.9267,
                           'Edxz': 27.9267,
                           'Edyz': 27.9267,
                           'Edxy': 27.9267,
                           'Edx2y2': 27.9267,
                           'Estar': 38.2661}}   
                                           

interactive_constans5 = {('C', 'C'): {'V_sssigma': V_alpha(h=-4.3882, r=1.42),
                                     'V_spsigma': V_alpha(h=5.4951, r=1.42),
                                     'V_sdsigma': V_alpha(h=-2.7655, r=1.42),
                                     'V_starsssigma':V_alpha(h=-2.3899, r=1.42),
                                     'V_starssigma': V_alpha(h=-2.6737, r=1.42),
                                     'V_starpsigma': V_alpha(h=5.1709, r=1.42),
                                     'V_stardsigma': V_alpha(h=-2.3034, r=1.42),
                                     'V_ppsigma': V_alpha(h=7.5480, r=1.42),
                                     'V_pppi': V_alpha(h=-2.6363, r=1.42),
                                     'V_pdsigma': V_alpha(h=-2.1621, r=1.42),
                                     'V_pdpi': V_alpha(h=3.9281, r=1.42),
                                     'V_ddsigma': V_alpha(h=-4.1813, r=1.42),
                                     'V_ddpi': V_alpha(h=4.9779, r=1.42),
                                     'V_ddd': V_alpha(h=-3.9884, r=1.42)}
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
                                          

configuration = {

'parametrization1':
{'lattice_type': 'rectagonal1_40_proc_of_defects',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1.42, # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 120,
 'horizontal_num_of_steps': 8,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies1,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 3600, #4910
 'interactive_constans': interactive_constans1,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': None,
 'number_of_defects': None},

 'parametrization2':
{'lattice_type': 'rectagonal2_40_proc_of_defects',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1.536329, # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 120,
 'horizontal_num_of_steps': 8,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies2,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 3600, #4910
 'interactive_constans': interactive_constans2,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': None,
 'number_of_defects': None},

 'parametrization3':
{'lattice_type': 'rectagonal3_40_proc_of_defects',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1.312, # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 120,
 'horizontal_num_of_steps': 8,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies3,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 3600, #4910
 'interactive_constans': interactive_constans3,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': None,
 'number_of_defects': None},

 'parametrization4':
{'lattice_type': 'rectagonal4_40_proc_of_defects',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1., # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 20,
 'horizontal_num_of_steps': 60,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies1,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 3600, #4910
 'interactive_constans': interactive_constans1,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': None,
 'number_of_defects': None},

 'parametrization5':
{'lattice_type': 'rectagonal5_40_proc_of_defects',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1.536329, # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 20,
 'horizontal_num_of_steps': 60,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies2,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 3600, #4910
 'interactive_constans': interactive_constans2,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': None,
 'number_of_defects': None},

 'parametrization6':
{'lattice_type': 'rectagonal6_40_proc_of_defects',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1.312, # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 20,
 'horizontal_num_of_steps': 60,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies3,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 4000, #4910
 'interactive_constans': interactive_constans3,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': None,
 'number_of_defects': None},

 'parametrization7':
{'lattice_type': 'rectagonal7_40_proc_of_defects',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1.42, # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 20,
 'horizontal_num_of_steps': 60,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies4,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 4000, #4910
 'interactive_constans': interactive_constans4,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': None,
 'number_of_defects': None},


 'parametrization8':
{'lattice_type': 'rectagonal8_40_proc_of_defects',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1.42, # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 120,
 'horizontal_num_of_steps': 8,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies4,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 3600, #4910
 'interactive_constans': interactive_constans4,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': None,
 'number_of_defects': None},


'parametrization9':
{'lattice_type': 'rectagonal9',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1.42, # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 20,
 'horizontal_num_of_steps': 60,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies5,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 4000, #4910
 'interactive_constans': interactive_constans5,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': None,
 'number_of_defects': None},


 'parametrization10':
{'lattice_type': 'rectagonal10',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1.42, # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 120,
 'horizontal_num_of_steps': 8,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies5,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 4000, #4910
 'interactive_constans': interactive_constans5,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': None,
 'number_of_defects': None},


 'parametrization11':
{'lattice_type': 'rectagonal11_40_proc_of_defects',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1.42, # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 20,
 'horizontal_num_of_steps': 60,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies5,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 4000, #4910
 'interactive_constans': interactive_constans5,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': True,
 'number_of_defects': None},


 'parametrization12':
{'lattice_type': 'rectagonal12_40_proc_of_defects',
 'gauss_sigma': 0.015,
 'start': -17,
 'stop': 17,  
 'step': 0.001,
 'distance': 1.42, # 1.536329 pomyśleć z tymi zmiennoprzecinkowymi 
 'lanczos_vectors':None,
 'vertical_num_of_steps': 120,
 'horizontal_num_of_steps': 8,
 'x_num_of_steps': 39,
 'sigma': 0.0001, 
 'magnitude': 'LM',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies5,
 'calculation_type': 'non spin',
 'number_of_eigenvalues': 4000, #4910
 'interactive_constans': interactive_constans5,
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',
 'ld':None,
 'lp': None,
 'number_of_friends':None,
 'defects': True,
 'number_of_defects': None}


} 


# TODO policzyc uklady z parametryzacją publikacje 
# TODO poiczyć dla dewastacji na brzegach heksagonalnej struktury i prostokątnej
#sys.argv[1]