import numpy as np 
import os 
import pandas as pd 
import glob
import natsort
from scipy.sparse.linalg import eigs, LinearOperator
import numpy as np
import scipy.sparse as sp

#matrix = scipy.sparse.load_npz('/home/przemek/Documents/Modeling/tight_binding/graphene_30_20/350_350.npz' )

path = '/home/przemek/Documents/Modeling/tight_binding/graphene_30_20/'
# TODO zadbać o odpowiednie łączenie macierzy plus zapamietywanie poprzedniego kroku
def get_batches(positions, dim):

    diffs = np.array([position - positions[i - 1] for i, position in enumerate(positions)][1:])
    split = np.where(diffs==dim)[0] + 1
    split = list(np.insert(split, 0, 0))
    split.append(None)
    indx = [(split[i], split[i+1]) for i in range(len(split) - 1)]
    batches = [positions[start_idx:end_idx] for start_idx, end_idx in indx]

    return batches

def load_matrices_positions(path):
    
    files = os.listdir(path)
    files = natsort.natsorted(files,reverse=False)
   
    matrix_files = []
    for num, f in enumerate(files):
    
        matrix_file = f.split('/', -1)[-1].split('.',-1)[0].split('.',-1)[0]
        matrix_file = matrix_file.split('_')
        
        matrix = np.load(path + f)
        matrix = sp.csr_matrix((matrix['data'], matrix['indices'] , matrix['indptr']), shape=matrix['shape'])
        #matrix_data = matrix['data']
        #matrix_indices = matrix['indices']
        #matrix_rows = np.array([idx for idx in matrix_indices]) + int(matrix_file[0])
        #matrix_columns = np.array([]) + int(matrix_file[1])
        #clprint(matrix_indices)
        #print(matrix.indices)
        print()
        print(matrix)
        
        
    return matrix_reshaped_indices
'''
def get_sorted_matrix_files(batches, path):
    
    files = os.listdir(path)
    files = natsort.natsorted(files,reverse=False)
    
    matrix_batches = [[path + next((s for s in files if str(position) in s), None)  for position in batch ] for batch in batches]
   

    return matrix_batches
    
def get_final_sparse_matrix(sorted_files):
    

    for sorted_file in sorted_files:
        horizontal_matrices = []
        for f in sorted_file: 
            print(f)
            matrix = np.load(f)
            matrix = sp.csr_matrix((matrix['data'], matrix['indices'] , matrix['indptr']), shape=matrix['shape'])
            print(matrix)

            #horizontal_matrices.append(matrix)
            
        #stacked_horizontally = sp.hstack(horizontal_matrices)
        #vertical_matrices.append(stacked_horizontally)
  

    return

'''
        
#zapytac prof co robic co bedzie lepiej utrzymywać :P bo mozna sparse matrices albo taki wektor , dodac gestośc stanów i to gotowe :p
print(path)
m = load_matrices_positions(path)
batches = get_batches(m, 10)
sorted_files = get_sorted_matrix_files(batches, path)

        #matrix = np.load(path + f)
        #matrix = sp.csr_matrix((matrix['data'], matrix['indices'] , matrix['indptr']),
         #            shape=matrix['shape']) 
'''
Xa = sp.csr_matrix((4, 4))
Xb = sp.csr_matrix((3, 6))


diff_n_rows = Xa.shape[0] - Xb.shape[0]
diff_n_columns = np.abs(Xa.shape[1] - Xb.shape[1])
print(diff_n_columns)

Xb_new_rows = sp.vstack((Xb, sp.csr_matrix((diff_n_rows, Xb.shape[1])))) 

X = sp.hstack((Xa, Xb_new_rows))

print('Xa',Xa.shape, 'Xb',Xb.shape, 'Xb_new_rows',Xb_new_rows.shape, 'X',X.shape)

l = [0, 0, 1, 1, 1, 2,2, 3,3, 4,4, 5, 5, 6,6, 7, 7, 8,8, 9,9, 10, 10, 10]

vstack = []
hstack = []

temp = l[0]
v = 0
for i in l:

    if temp == i:
        print('vstack', i)
        vstack.append(i)
        v = v + i 

    
    
    #print('v', v)
    
    

    if temp != i:
        print('h', i, 'v',v)
        hstack.append(i)
        temp = i

    
print('hstack',hstack, 'vstack',vstack)

    temp_matrix_file = matrix_files[0]
    matrix = np.load(path + files[0])
    vstack = sp.csr_matrix((matrix['data'], matrix['indices'] , matrix['indptr']),
                     shape=matrix['shape']) 
    
    for matrix_file, f in zip(matrix_files[1:], files[1:]):
    
        matrix = np.load(path + f)
        matrix = sp.csr_matrix((matrix['data'], matrix['indices'] , matrix['indptr']),
                     shape=matrix['shape'])  

        if temp_matrix_file == matrix_file:
            
            
            if vstack.shape[1] != matrix.shape[1]:
            
                diff_n_rows = vstack.shape[1] - matrix.shape[1]
                matrix_new = sp.hstack((matrix, sp.csr_matrix((diff_n_rows, matrix.shape[0]))))
                print(matrix_file, vstack.shape, matrix.shape, matrix_new.shape)
            
            else:
                matrix_new = matrix 
            
            vstack = sp.vstack([vstack, matrix_new])
                
        elif temp_matrix_file != matrix_file:
            
            if vstack.shape[0] != matrix.shape[0]:
            
                diff_n_rows = vstack.shape[0] - matrix.shape[0]
                matrix_new = sp.vstack((matrix, sp.csr_matrix((diff_n_rows, matrix.shape[1]))))
            
            else:
                matrix_new = matrix
                
            vstack = sp.hstack([vstack, matrix_new])
            temp_matrix_file = matrix_file
        '''
