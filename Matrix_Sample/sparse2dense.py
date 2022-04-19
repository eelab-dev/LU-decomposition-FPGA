import numpy as np
from scipy.sparse import csc_matrix
from scipy.linalg import lu
import scipy.io as io
import os

os.chdir(os.path.dirname(__file__))

"""
indptr = np.array([0, 2, 4, 6, 7, 8])
indices = np.array([0, 3, 1, 2, 2, 3, 3, 4])
data = np.array([1, 1, 1, 7 / 9, 1, 0, 1, 1])
L = csc_matrix((data, indices, indptr), shape=(5, 5)).toarray()
print(L)

indptr = np.array([0, 1, 2, 4, 7, 8])
indices = np.array([0, 1, 0, 2, 1, 2, 3, 4])
data = np.array([1, 1 / 7, 0.5, 5 / 18, 1, 2 / 9, 1 / 6, 1])
U = csc_matrix((data, indices, indptr), shape=(5, 5)).toarray()
print(U)

indptr = np.array([0, 0, 0, 0, 0, 1])
indices = np.array([2])
data = np.array([1 / 9])
F = csc_matrix((data, indices, indptr), shape=(5, 5)).toarray()
print(F)
"""

# data = np.array(io.loadmat("matlab.mat")["A"])
data = np.array([[-7, 0, 0, 7, -8], [-10, 0, -7, 0, 0], [-1, -8, 1, 10, 0], [0, 9, -10, 0, 6], [2, 0, 6, -5, -9]])
print(data)
A = csc_matrix(data)
print(repr(A.indptr), repr(A.indices), repr(A.data), sep='\n')
print(repr(data.flatten()))
P, L, U = lu(data)
print(P, L, U, sep='\n')
