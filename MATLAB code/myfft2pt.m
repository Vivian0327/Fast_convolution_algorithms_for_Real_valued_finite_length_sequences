%% myfft2pt

% 7/17/2023 FAMU-FSU colledge of engineering
% This file is coded by Rajesh Thomas and Weiwei Wang (ww20br@fsu.edu).
% Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu).
% Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

% myfft2pt: This function works by recursively prioritizing 2-point DFTs.
%   If the length is odd and is less than gpfact, the DFT is computed using 
%   Matveev's eigenvectors.

% Syntax:
%   [y a] = myfft2pt(x)

% Description:
%   [y a] = myfft2pt(x) returns the DFT of the column vector x. 

% Input Arguments:
%   x: Input signal in a column vector.

% Output Arguments:
%   y: the normalized DFT of x.
%   a: is a two elements column vector and it's a global argument used to 
%      store the number of operations.a(1) is the number of
%      multiplications, a(2) is the number of addtions.



function [y a]= myfft2pt(x)
global a; % a stores the number of operations.
if isempty(a)
    a = [0; 0];
end

N = size(x,1);
gpfact = 30;  % greatest possible factor. reduce gpfact for more accuracy.
chkfact = 1:gpfact;
fact = chkfact(rem(N,chkfact)==0); % factors of N <= gpfact

y = zeros(N, 1, 'like', 1j);
if N == 2
    y(1) = x(1) + x(2);
    y(2) = x(1) - x(2);
    y = y/sqrt(2);
    if isreal(x)
        a = a+[0; 2];
    else
        a = a+[0; 4];
    end
elseif fact(end) == 1 || (rem(N,2)==1 && N <= gpfact)
    % i.e., if N is prime for factors <= gfact or if
    % N is odd and N <= gpfact
    y = dftf(x);
    a = a + numops(N,isreal(x));
else
    m = N/fact(2); % m is length of DFT block
    if fact(2) ~= 2
        m = N/fact(end);
        if m <= gpfact
            m = fact(end);
        end
    end
    z = zeros(N, 1, 'like', 1j);
    w = exp(-1j*2*pi/N*(0:m-1)');
    for k = 1:N/m  % loop for m-point DFT blocks
        z((k-1)*m+1:k*m) = myfft2pt(x(k:N/m:N)).*(w.^(k-1));
    end
    if N ~= 4
        a = a + (N/m-1)*[4*(m-1); 2*(m-1)];
    end
    
    for k = 1:m  % final N/m-point DFT butterflies
        y(k:m:N) = myfft2pt(z(k:m:N));
    end
end
end
