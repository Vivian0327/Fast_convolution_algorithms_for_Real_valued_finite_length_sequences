'''
primeFactors

7/17/2023 FAMU-FSU colledge of engineering
Coded by Weiwei Wang(ww20br@fsu.edu). 
Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu). 
Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

primeFactors: 
    A function to find all prime factors of a given number n.
parameters:
    n: a scalar
       Input parameter n
Returns:
    prime_factors: array_like
                   the prime factors of input n

This function is called by "myfftpfa".
'''

import math
 
def primeFactors(n):
     
    i = 2
    prime_factors = []
    while i*i <= n:
        if n%i == 0:
            prime_factors.append(i)
            n //= i
        else:
            i += 1
    
    if n>1:
        prime_factors.append(n)
    
    return prime_factors
 