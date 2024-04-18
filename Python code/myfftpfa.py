'''
myfftpfa

7/17/2023 FAMU-FSU colledge of engineering
Coded by Weiwei Wang(ww20br@fsu.edu). 
Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu). 
Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

myfftpfa: ufunc
          calculates the fft using Prime Factor Algorithm (PFA). Doing the fft need to call function "myfft2pt".

Description:
   myfftpfa(x) returns the fft of the input signal x. The output y is the fft of x.

Parameters:
   x: array_like
      Input signal in a array.

Returns:
   y: array_like
      the fft of x.

Reference:
   https://urldefense.com/v3/__https://en.wikipedia.org/wiki/Prime-factor_FFT_algorithm__;!!Epnw_ITfSMW4!6l3ZUXxXrAsmc9JaL9TqDfRq_R3O0QLXddTCWJNQG8r4D_6IDn5DHaHupEDf8kahywg$ 

'''

import numpy as np
from primefactors import primeFactors
from myfft2pt import myfft2pt

# Python program for the extended Euclidean algorithm
def gcd(a, b):
    if b == 0:
        gcd, s, t = a, 1, 0
        return (gcd, s, t)
    else:
        s2, t2, s1, t1 = 1, 0, 0, 1
    while b > 0:
        q= a // b
        r, s, t = (a - b * q),(s2 - q * s1),( t2 - q * t1)
        a,b,s2,t2,s1,t1=b,r,s1,t1,s,t
        
    gcd,s,t=a,s2,t2
    
    return (gcd,s,t)


def myfftpfa(x):
    
    if not isinstance(x, np.ndarray) or x.ndim != 2 or x.shape[1] != 1:
        raise ValueError('Input should be a column vector.')
        
    gpfact = 25             # refer myfft2pt() for more info on gpfact
    N = len(x)
    fact = np.array(np.unique(primeFactors(N)))  # prime factors of N
    ufac = np.unique(fact)  # unique prime factors of N
    nfac = len(ufac)        # number of unique prime factors
    
    y = np.zeros((N, 1), dtype=complex)
    
    if nfac == 1:           # length of input is a power of a prime number
        y = myfft2pt(x)
    else:
        z = np.zeros((N, 1), dtype=complex)
        powfac = np.zeros(nfac, dtype=int)  # powers of unique prime factors of N
        for i in range(nfac):
            powfac[i] = np.sum(fact == ufac[i])
        
        N1 = ufac[0] ** powfac[0]           # N1 * N2 = N; N1, N2 are co-prime
        N2 = N // N1
        
        if ufac[0] == 3 and powfac[0] >= 2:  # 9 is a special case...
            # power of a prime number other than 2 (i.e., 3^2) is less than gpfact.
            temp = ufac[ufac < gpfact]
            if temp[-1] < 9:                 # if the largest prime factor < 3^2
                N2 = ufac[0] ** powfac[0]
                N1 = N // N2
        
        N2idx = np.arange(N2) * N1
        for m in range(N1):                 # loop to find N2-point DFTs
            idx = np.mod(N2idx + m * N2, N) + 1
            z[m * N2:(m + 1) * N2] = np.c_[myfftpfa(x[idx - 1])]
        
        # this works only if N1 and N2 are co-prime
        _, N1inv, N2inv = gcd(N1, N2)      # N1 * N1inv = 1 (mod N2); N2 * N2inv = 1 (mod N1)
        
        K1idx = np.arange(N1) * N2inv * N2
        for m in range(N2):                # loop to find N1-point DFTs
            idx = np.mod(K1idx + m * N1inv * N1, N) + 1
            y[idx - 1] = np.c_[myfftpfa(z[m::N2])]
            
    return y
