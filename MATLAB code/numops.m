%% numops

% 7/17/2023 FAMU-FSU colledge of engineering
% This file is Coded by Rajesh Thomas modified by Weiwei Wang (ww20br@fsu.edu).
% Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu).
% Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

% numops: Calculate the number of multiplications and number of additions 
%   for calculating the DFT using Matveev's DFT eigenvectors.

% Syntax:
%   ops = numops(n, type) 

% Description:
%   ops = numops(n, type) return the number of multiplications and the
%   number of additions for calculating the n-point DFT.

% Input Arguments:
%   n = length of DFT
%   type = type of DFT input. 'r' = real, 'c' = complex

% Output Arguments:
%   ops(1) = number of multiplications
%   ops(2) = number of additions

function ops = numops(n, type)

if ~isscalar(n)
    error('Length of DFT must be a scalar.');
end

m = floor([n+4, n+2, n-1, n+1]/4);

mp1 = (m(1)-1)*(floor(n/2)-(m(1)-2)/2) + 2;
mn1 = (m(2)-1)*(floor(n/2)-(m(2)-2)/2) + 2;
mpj = m(3)*(floor((n-1)/2)-(m(3)-1)/2);
mnj = m(4)*(floor((n-1)/2)-(m(4)-1)/2);
mul = 2*(mp1+mn1+mpj+mnj);

ap1x = m(1)*n - m(1)*(m(1)-1) - 1;
an1x = m(2)*n - m(2)*(m(2)-1) - 1;
apjx = 2*m(3)*floor((n-1)/2) - m(3)^2;
anjx = 2*m(4)*floor((n-1)/2) - m(4)^2;

apn1w = m(2)^2 + (m(1)+m(2)-1)*(floor(n/2) + 1 - m(2));
apnjw = m(3)^2 + (m(3)+m(4)-1)*(floor((n-1)/2) - m(3));

add = (ap1x+an1x+apjx+anjx) + (apn1w+apnjw);

 % if the DFT input is complex
if type == 'c' || type == 0   
    mul = 2*mul;
    add = 2*add + 2*n;
end

ops = [mul; add];

end