'''
numops

7/17/2023 FAMU-FSU colledge of engineering
Coded by Weiwei Wang(ww20br@fsu.edu). 
Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu). 
Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

numops: 
    Calculate the number of multiplications and number of additions for calculating the DFT using Matveev's DFT eigenvectors.

Description:
   numops(n, type) return the number of multiplications and the
   number of additions for calculating the n-point DFT.

Parameters:
   n = length of DFT
   type = type of DFT input. 'r' = real, 'c' = complex

Returns:
   ops: array_like
       ops(1) = number of multiplications
       ops(2) = number of additions

'''

def numops(n, type):
    
    #if not isinstance(n, int):
    #    raise ValueError('Length of DFT must be an integer.')
    
    m = [int((n + 4) / 4), int((n + 2) / 4), int((n - 1) / 4), int((n + 1) / 4)]
    
    mp1 = (m[0] - 1) * (int(n / 2) - int((m[0] - 2) / 2)) + 2
    mn1 = (m[1] - 1) * (int(n / 2) - int((m[1] - 2) / 2)) + 2
    mpj = m[2] * (int((n - 1) / 2) - int((m[2] - 1) / 2))
    mnj = m[3] * (int((n - 1) / 2) - int((m[3] - 1) / 2))
    mul = 2 * (mp1 + mn1 + mpj + mnj)
    
    ap1x = m[0] * n - m[0] * (m[0] - 1) - 1
    an1x = m[1] * n - m[1] * (m[1] - 1) - 1
    apjx = 2 * m[2] * int((n - 1) / 2) - m[2] ** 2
    anjx = 2 * m[3] * int((n - 1) / 2) - m[3] ** 2
    
    apn1w = m[1] ** 2 + (m[0] + m[1] - 1) * (int(n / 2) + 1 - m[1])
    apnjw = m[2] ** 2 + (m[2] + m[3] - 1) * (int((n - 1) / 2) - m[2])
    
    add = (ap1x + an1x + apjx + anjx) + (apn1w + apnjw)
    
    if type == 'c' or type == 0:
        mul = 2 * mul
        add = 2 * add + 2 * n
    
    ops = [mul, add]
    
    return ops
