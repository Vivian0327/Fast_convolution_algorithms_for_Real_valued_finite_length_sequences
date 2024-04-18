'''
dfteig

7/17/2023 FAMU-FSU colledge of engineering
Coded by Weiwei Wang(ww20br@fsu.edu). 
Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu). 
Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

dfteig：
    Caluculate the eigenvectors and the multiplicity of the eigenvalues.

Docstring:
   V = dfteig(k) returns a square matrix V whose columns contain the 
   real and orthonormal eigenvectors of the k-point DFT. These eigenvectors
   are computed using the Matveev's method.

   V, emult = dfteig(k) also returns a 4 length row vector emult,
   which contains the multiplicity of the DFT eigenvalues. The
   multiplicity corresponds to the eigenvalues {1,-1,j,-j} in the
   same order. Thus, emult(1) is the eigenvalue multiplicity for the
   eigenvalue 1, emult(2) is the multiplicity for -1, and so on.

   dfteig(k,'inv') returns the orthonormal eigenvectors and/or the
   eigenvalue multiplicity for the inverse Fourier transform. The
   eigenvalue multiplicity is in the same order as the forward
   transform, i.e., {1,-1,j,-j}.

Parameters:
   k: the length of the input signal or the length of the DFT
   'inv': means the inverse Fourier transform

Returns:
   V: is a square(k*k) matrix whose columns are the real and orthonormal
      eigenvectors
   emult: is a row matrix with the length 4. The elements of emult are the
      multiplicity corresponds to the eigenvalues {1,-1,j,-j} in the same
      order.
'''

import numpy as np
import math

def dftmtx(N):
    #return sp.fft(sp.eye(N))
    return np.fft.fft(np.eye(N))

def gramDet(p, evenodd):
    k, mult = p.shape
    
    if mult == 0:
        return np.zeros((k, 0))

    v = np.column_stack((p[:, 0]/np.linalg.norm(p[:, 0]), np.zeros((k,mult-1))))
    
    for i in range(2, mult+1):
        sign = (-1)**(i-1)
        ptem = p[evenodd:(i+evenodd), :(i-1)]
        
        for j in range(1, i+1):
            pptem = ptem[np.concatenate((np.arange(j-1), np.arange(j, i))), :]
            cofactor = np.linalg.det(pptem)
            
            v[:, i-1] = v[:, i-1] + sign * cofactor * p[:, j-1]
            sign = sign * (-1)
        
        v[:, i-1] = v[:, i-1]/np.linalg.norm(v[:,i-1])
        
    return v

def dfteig(k, *args):

    if k <= 0 or k != int(k):
        raise ValueError('k should be an integer greater than 0.')

    # F is the Fourier transform matrix
    F = dftmtx(k)/math.sqrt(k);   

    # emult is eigenvalue multiplicity in the order {1, -1, j, -j}.
    emult = np.floor(np.multiply([k+4, k+2, k-1, k+1],1/4))
    
    # = F^2
    F2 = np.fliplr(np.eye(k-1))
    F2_2 = np.array(np.zeros((1,k-1))) 
    F2 = np.insert(F2,0,F2_2,0)
    F2_2 = np.array(np.zeros((1,k))) 
    F2 = np.insert(F2,0,F2_2,1)
    F2[0,0]=1
    F2[0,k-1]=0

    # conjugate then transpose
    F3 = np.conj(F).T
    
    # Projection matrices
    pp1 = (F3  + F2 + F + np.eye(k))/4
    pn1 = (-F3 + F2 - F + np.eye(k))/4
    ppj = (1j*F3  - F2 - 1j*F + np.eye(k))/4
    pnj = (-1j*F3 - F2 + 1j*F + np.eye(k))/4
    
    emult1 = int(emult[0])
    emult2 = int(emult[1])
    emult3 = int(emult[2]+1)
    emult4 = int(emult[3]+1)
    
    pp1 = np.real(pp1[:,0:emult1])
    pn1 = np.real(pn1[:,0:emult2])
    ppj = np.real(ppj[:,1:emult3])
    pnj = np.real(pnj[:,1:emult4])
    
    # even = 0 (i.e., even projection matrix)
    vp1 = gramDet(pp1, 0) 
    vn1 = gramDet(pn1, 0)
    
    # odd = 1 (i.e., odd projection matrix)
    vpj = gramDet(ppj, 1) 
    vnj = gramDet(pnj, 1)
    
    v = np.hstack((vp1, vn1, vpj, vnj)) # orthonormal eigenvectors

    if len(args) == 1 and args[0] == 'inv':
        emult[2:4] = np.floor([k + 1, k - 1] / 4)
        v = np.hstack((vp1, vn1, vnj, vpj))
    elif len(args) > 1:
        raise ValueError('There should be a maximum of only 2 input arguments.')

    if np.sum(np.array(v.shape) == k) != 2:
        raise ValueError('Wrong code.')
        
    return v, emult