'''
myfft2pt

7/17/2023 FAMU-FSU colledge of engineering
Coded by Weiwei Wang(ww20br@fsu.edu). 
Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu). 
Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

myfft2pt:
    This function works by recursively prioritizing 2-point DFTs.If the length is odd and is less than gpfact,the DFT is computed using Matveev's eigenvectors.
    This function is called by "myfftpfa" to do the fft and can calculate the number of operations. 

Description:
    myfft2pt(x) returns the DFT of the input signal x. 

Parameters:
   x: array_like
      Input signal in a array.

Returns:
   y: array_like
      the DFT of x.
   a: global variable
      is a two elements column vector and it's a global argument used to 
      store the number of operations.a(1) is the number of
      multiplications, a(2) is the number of addtions.
 
'''

import numpy as np

from dftf import dftf
from dfteig import dfteig
from numops import numops

import globals

def myfft2pt(x):

    N = x.shape[0]
    gpfact = 30  # greatest possible factor. reduce gpfact for more accuracy.
    chkfact = np.arange(1, gpfact + 1)
    fact = chkfact[N % chkfact == 0]  # factors of N <= gpfact
    
    y = np.zeros((N, 1), dtype=complex)
    
    if N == 2:
        y[0] = x[0] + x[1]
        y[1] = x[0] - x[1]
        y /= np.sqrt(2)
        
        if np.isrealobj(x):
            globals.a += np.array([0, 2], dtype=int)
        else:
            globals.a += np.array([0, 4], dtype=int)
    elif fact[-1] == 1 or (N % 2 == 1 and N <= gpfact):
        # i.e., if N is prime for factors <= gfact or if
        # N is odd and N <= gpfact
        y = dftf(x)
        globals.a += numops(N, np.isrealobj(x))
    else:
        m = N // fact[1]  # m is length of DFT block
        
        if fact[1] != 2:
            m = N // fact[-1]
            if m <= gpfact:
                m = fact[-1]
        
        z = np.zeros((N, 1), dtype=complex)
        w = np.c_[np.exp(-1j * 2 * np.pi / N * np.arange(m))]
        
        # loop for m-point DFT blocks
        for k in range(1, N // m + 1):         
            z[(k - 1) * m:k * m] = myfft2pt(x[k - 1:N:N // m]) * (w ** (k - 1))
        
        if N != 4:
            globals.a += (N // m - 1) * np.array([4 * (m - 1), 2 * (m - 1)], dtype=int)
        
        # final N/m-point DFT butterflies
        for k in range(1, m + 1):  
            y[k - 1:m:N] = myfft2pt(z[k - 1:m:N])
    
    return np.c_[y]
