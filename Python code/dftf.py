'''
dftf

7/17/2023 FAMU-FSU colledge of engineering
Coded by Weiwei Wang(ww20br@fsu.edu). 
Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu). 
Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

dftf: 
   calculates the discrete Fourier transform (DFT) using its real and orthonormal eigenvectors.

Docstring:
   dftf(x) returns the DFT of the column vector x. The output y is normalized such that the DFT is unitary, i.e., y = fft(x)/sqrt(length(x)).

Parameters:
   x: array_like
      Input signal x.

returns:
   y: array_like
      the normalized DFT of x.

   The eigenvectors are calculated via the dfteig() function which uses Matveev's method to get the orthonormal eigenvectors.

Call signature:
    y = dftf(x)
'''

import torch
import numpy as np
import scipy as sp

from dfteig import dfteig

def dftf(x):
    
    #if not np.ndim(x) == 1:
    #    raise ValueError('Input should be a column vector.')
        
    n = len(x)
    v, m = dfteig(n)
   
    idx = np.cumsum(m)
    idx1 = int(idx[1])
    
    # combiner matrix
    u = np.zeros((n, 2))
    u[:idx1, 0] = np.concatenate((np.ones(int(m[0])), -np.ones(int(m[1]))))
    u[idx1:, 1] = np.concatenate((np.ones(int(m[2])), -np.ones(int(m[3]))))
    
    xtilde = np.diagflat( v.T @ x) # output after filtering and downsampling
    v_scale = v @ xtilde
    y = v_scale @ u                     
    y = y @ np.array([1, 1j])  
    y = np.sqrt(n) * y

    return y