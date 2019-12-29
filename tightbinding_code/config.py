import numpy as np
def V_alpha(h, r=1.42, type_structure='graphene'):
  n = 2.
  n_c = 6.5
  r_o = 2.23
  r_c = 2.18
  if type_structure == 'graphene':
    s = np.exp( -n * (r / r_c)**(n_c) + n * (r_o / r_c)**(n_c) ) * (r_o / r)**n 
  else:
    s = (r_o / r)**n 
  h_alpha = h * s
  return h_alpha

diagonal_energies_p = {'C': {'Es': -50,
                           'Epx': -50,
                           'Epy': -50,
                           'Epz': 0,
                           'Edz2': 4000.,
                           'Edxz': 4000.,
                           'Edyz': 4000.,
                           'Edxy': 4000.,
                           'Edx2y2': 4000.,
                           'Estar': 4000.}}   

                              
interactive_constans_p = {('C', 'C'): {'V_sssigma': 0,
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



diagonal_energies_sp = {'C': {'Es': -2.99,
                           'Epx': 3.71,
                           'Epy': 3.71,
                           'Epz': 3.71,
                           'Edz2': 4000.,
                           'Edxz': 4000.,
                           'Edyz': 4000.,
                           'Edxy': 4000.,
                           'Edx2y2': 4000.,
                           'Estar': 4000.}}   
                                           

interactive_constans_sp = {('C', 'C'): {'V_sssigma': V_alpha(h=-5., r=1.42),
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




diagonal_energies_spds = {'C': {'Es': -1.0458,
                           'Epx': 7.0850,
                           'Epy': 7.0850,
                           'Epz': 7.0850,
                           'Edz2': 27.9267,
                           'Edxz': 27.9267,
                           'Edyz': 27.9267,
                           'Edxy': 27.9267,
                           'Edx2y2': 27.9267,
                           'Estar': 38.2661}}   
                                           

interactive_constans_spds = {('C', 'C'): {'V_sssigma': V_alpha(h=-4.3882, r=1.42),
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


diagonal_energies_silicen = {'C': {'Es': -2.0196,
                           'Epx': 4.5448,
                           'Epy': 4.5448,
                           'Epz': 4.5448,
                           'Edz2': 14.1836,
                           'Edxz': 14.1836,
                           'Edyz': 14.1836,
                           'Edxy': 14.1836,
                           'Edx2y2': 14.1836,
                           'Estar': 19.6748}}   
                                           

interactive_constans_silicen = {('C', 'C'): {'V_sssigma': V_alpha(h=-1.9413, r=1.42, type_structure='silicen'),
                                     'V_spsigma': V_alpha(h=2.7836, r=1.42, type_structure='silicen'),
                                     'V_sdsigma': V_alpha(h=-2.7998, r=1.42, type_structure='silicen'),
                                     'V_starsssigma':V_alpha(h=-1.6933, r=1.42, type_structure='silicen'),
                                     'V_starssigma': V_alpha(h=-3.3081, r=1.42, type_structure='silicen'),
                                     'V_starpsigma': V_alpha(h=2.8428, r=1.42, type_structure='silicen'),
                                     'V_stardsigma': V_alpha(h=-0.7003, r=1.42, type_structure='silicen'),
                                     'V_ppsigma': V_alpha(h=4.1068, r=1.42, type_structure='silicen'),
                                     'V_pppi': V_alpha(h=-1.5934, r=1.42, type_structure='silicen'),
                                     'V_pdsigma': V_alpha(h=-2.1073, r=1.42, type_structure='silicen'),
                                     'V_pdpi': V_alpha(h=1.9977, r=1.42, type_structure='silicen'),
                                     'V_ddsigma': V_alpha(h=-1.2327, r=1.42, type_structure='silicen'),
                                     'V_ddpi': V_alpha(h=2.5145, r=1.42, type_structure='silicen'),
                                     'V_ddd': V_alpha(h=-2.4734, r=1.42, type_structure='silicen')}
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

settings = { 
 'gauss_sigma': 0.015,
 'start': -10,
 'stop': 20,  
 'step': 0.001,
 'title': 'hexagonal_ZIGZAG_silicen',
 'saving_directory': '/home/przemek/Documents/Modeling/tight_binding/results_diploma',}
                                          

configuration = {

'parametrization':
{
 
 'ld':None,
 'lp': None,
 'magnitude': 'LM',
 'distance': 1.42,
 'x_num_of_steps': 3,
 'lanczos_vectors':None,
 'number_of_friends':None,
 'number_of_eigenvalues': 1,
 'calculation_type': 'non spin',
 'neighbour_calculation_method': 'distance',
 'diagonal_energies':diagonal_energies_silicen,
 'interactive_constans': interactive_constans_silicen,
 'fermi_level': 0
  }} 
