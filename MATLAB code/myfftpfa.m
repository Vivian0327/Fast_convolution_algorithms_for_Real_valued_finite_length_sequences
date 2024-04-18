%% myfftpfa

% 7/17/2023 FAMU-FSU colledge of engineering
% This file is Coded by Rajesh Thomas modified by Weiwei Wang (ww20br@fsu.edu).
% Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu).
% Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

% myfftpfa: calculates the fft using Prime Factor Algorithm (PFA).

% Syntax:
%   y = myfftpfa(x)

% Description:
%   y = myfftpfa(x) returns the fft of the column vector x. The output y is
%       the result of fft.

% Input Arguments:
%   x: Input signal in a column vector.

% Output Arguments:
%   y: the fft of x.

% Reference:
%   https://urldefense.com/v3/__https://en.wikipedia.org/wiki/Prime-factor_FFT_algorithm__;!!Epnw_ITfSMW4!6l3ZUXxXrAsmc9JaL9TqDfRq_R3O0QLXddTCWJNQG8r4D_6IDn5DHaHupEDf8kahywg$ 

function y = myfftpfa(x)
if ~iscolumn(x)
    error('Input should be a column vector.');
end

gpfact = 25;         % refer myfft2pt() for more info on gpfac
N = length(x);
fact = factor(N);    % prime factors of N
ufac = unique(fact); % unique prime factors of N
nfac = length(ufac); % number of unique prime factors

y = zeros(N,1,'like',1j);
if nfac == 1         % length of input is a power of a prime number
    y = myfft2pt(x);
else
    z = zeros(N,1,'like',1j);
    powfac = zeros(1,nfac); % powers of unique prime factors of N
    for i = 1:nfac
        powfac(i) = sum(fact==ufac(i));
    end
    N1 = ufac(1)^powfac(1); % N1*N2=N; N1, N2 are co-prime
    N2 = N/N1;
    if ufac(1) == 3 && powfac(1) >= 2 % 9 is a special case...
                % power of a prime number other than 2 (ie, 3^2) is less than gpfact.
        temp = ufac(ufac < gpfact);
        if temp(end) < 9  % if the largest prime factor < 3^2
            N2 = ufac(1)^powfac(1);
            N1 = N/N2;
        end
    end    

    N2idx = (0:N2-1)*N1;
    for m = 0:N1-1 % loop to find N2-point DFTs
        idx = mod(N2idx + m*N2, N) + 1;
        z(m*N2+1:(m+1)*N2) = myfftpfa(x(idx));
    end
    
    [~, N1inv, N2inv] = gcd(N1, N2); % N1*N1inv = 1 (mod N2); N2*N2inv = 1 (mod N1)
                            % this works only if N1 and N2 are co-prime
    
    K1idx = (0:N1-1)*N2inv*N2;
    for m = 0:N2-1 % loop to find N1-point DFTs
        idx = mod(K1idx + m*N1inv*N1, N) + 1;
        y(idx) = myfftpfa(z(m+1:N2:N));
    end
end