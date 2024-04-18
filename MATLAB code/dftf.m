%% dftf: 
%   calculates the discrete Fourier transform (DFT) using its real and
%   orthonormal eigenvectors.

% 7/17/2023 FAMU-FSU colledge of engineering
% This file is coded by Rajesh Thomas and modified by Weiwei Wang (ww20br@fsu.edu).
% Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu).
% Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

% Syntax:
%   y = dftf(x)

% Description:
%   y = dftf(x) returns the DFT of the column vector x. The output y is
%   normalized such that the DFT is unitary, i.e.,
%   y = fft(x)/sqrt(length(x)).

% Input Arguments:
%   x: Input signal in a column vector.

% Output Arguments:
%   y: the normalized DFT of x.

%   The eigenvectors are calculated via the dfteig() function which uses
%   Matveev's method to get the orthonormal eigenvectors.

% Example:
% x = [-2; 0; 3; 1; 1; 3; 2];
% y = dftf(x)
% output:
% y =
%   8.0000 + 0.0000i
%  -3.8901 + 1.5637i
%  -6.6039 + 1.9499i
%  -0.5060 + 0.8678i
%  -0.5060 - 0.8678i
%  -6.6039 - 1.9499i
%  -3.8901 - 1.5637i

function y = dftf(x)

if ~iscolumn(x)
    error('Input should be a column vector.');
end
n = size(x,1);
[v, m] = dfteig(n);

idx = cumsum(m);

% combiner matrix
u = zeros(n,2); 
u(1:idx(2),1) = [ones(m(1),1); -ones(m(2),1)];
u(idx(2)+1:end,2) = [ones(m(3),1); -ones(m(4),1)];

xtilde = diag(v.'*x); % output after filtering and downsampling
v_scale = v*xtilde;
y = v_scale*u;
y = y*[1 1j].';
y = sqrt(n)*y;

end