%% dfteig

% 7/17/2023 FAMU-FSU colledge of engineering
% This file is proveded by Rajesh Thomas.
% Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu).
% Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

%   Syntax:
%      V = dfteig(k) 
%      [V, emult] = dfteig(k)
%      V = dfteig(k,'inv')
%      [V, emult] = dfteig(k,'inv')

%   Description:
%      V = dfteig(k) returns a square matrix V whose columns contain the 
%   real and orthonormal eigenvectors of the k-point DFT. These eigenvectors
%   are computed using the Matveev's method.
%
%      [V, emult] = dfteig(k) also returns a 4 length row vector emult,
%   which contains the multiplicity of the DFT eigenvalues. The
%   multiplicity corresponds to the eigenvalues {1,-1,j,-j} in the
%   same order. Thus, emult(1) is the eigenvalue multiplicity for the
%   eigenvalue 1, emult(2) is the multiplicity for -1, and so on.
%
%       dfteig(k,'inv') returns the orthonormal eigenvectors and/or the
%   eigenvalue multiplicity for the inverse Fourier transform. The
%   eigenvalue multiplicity is in the same order as the forward
%   transform, i.e., {1,-1,j,-j}.

% Input Arguments:
%   k: the length of the input signal or the length of the DFT
%   'inv': means the inverse Fourier transform

% Output Arguments:
%   V: is a square(k*k) matrix whose columns are the real and orthonormal
%      eigenvectors
%   emult: is a row vextor with the length 4. The elements of emult are the
%      multiplicity corresponds to the eigenvalues {1,-1,j,-j} in the same
%      order.

function [v, emult] = dfteig(k, varargin)

if k <= 0 || floor(k) ~= k
    error('k should be an integer greater than 0.');
end

F = dftmtx(k)/sqrt(k);   % F is the Fourier transform matrix
% F = symdhtmtx(k,k);

emult = floor([k+4, k+2, k-1, k+1]/4);
% emult is eigenvalue multiplicity in the order {1, -1, j, -j}.

F2 = eye(k); F2 = [F2(:,1) fliplr(F2(:,2:k))]; % = F^2
F3 = F'; % = F^3

pp1 = (F3  + F2 + F + eye(k))/4; % Projection matrices
pn1 = (-F3 + F2 - F + eye(k))/4;
ppj = (1j*F3  - F2 - 1j*F + eye(k))/4;
pnj = (-1j*F3 - F2 + 1j*F + eye(k))/4;

pp1 = real(pp1(:,1:emult(1)));
pn1 = real(pn1(:,1:emult(2)));
ppj = real(ppj(:,2:emult(3)+1));
pnj = real(pnj(:,2:emult(4)+1));

vp1 = gramDet(pp1, 0); % even = 0 (i.e., even projection matrix)
vn1 = gramDet(pn1, 0);
vpj = gramDet(ppj, 1); % odd = 1 (i.e., odd projection matrix)
vnj = gramDet(pnj, 1);

v = [vp1 vn1 vpj vnj]; % orthonormal eigenvectors

if nargin == 2 && isequal(varargin{1}, 'inv')
    emult(3:4) = floor([k+1, k-1]/4);
    v = [vp1 vn1 vnj vpj];
elseif nargin > 2
    error('There should be a maximum of only 2 input arguements.');
end

if sum(size(v) == k) ~= 2
    error('Wrong code.')
end

end

% function to calculate Grammian determinants.
function [v] = gramDet(p, evenodd)

[k, mult] = size(p);
if mult == 0
    v = zeros(k,0);
    return;
end

v = [p(:,1)/norm(p(:,1)) zeros(k, mult-1)];
for i = 2:mult
    sign = (-1)^(i-1);
    ptem = p((1:i)+evenodd,(1:i-1));
    for j = 1:i
        cofactor = det(ptem([1:j-1, j+1:i], :));
        v(:,i) = v(:,i) + sign*cofactor*p(:,j);
        sign = sign*(-1);
    end
    v(:,i) = v(:,i)/norm(v(:,i));
end
end