%% Real_conv_Computation_Fig2and3

% 7/17/2023 FAMU-FSU colledge of engineering
% Programmed by Weiwei Wang (ww20br@fsu.edu).
% Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu).
% Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

% This code is used to do the computation complexity comparasion for 
% Real-valued input signals.
% The length of sequence h is fixed.
% The length of sequence x is varied.

% 1. We need to call the function "opnums" to get the number of operations 
%    of the algorithm RV_DFT .
% 2. We need to call the function "myfftpfa->myfft2pt->dftf->dfteig" to get
%    the number of operations of the RV_FFT.

%% 
close all;
clear;

Lx = [26 45 65 102 115 130 189 240 256 500 728 1013 1150]; Lh = 12; J = 13;
%  Lx = [21 37 51 96 117 181 251 294 387 493 581 894 1002 1131]; Lh = 20; J = 21; 

clear global a;
global a;

% Prelocating the matrixes for speed
N = length(Lx);
RM_ORRVconv = zeros(1,N); RA_ORRVconv = zeros(1,N);
RM_RVconv = zeros(1,N);   RA_RVconv = zeros(1,N);
RM_DHTconv = zeros(1,N);  RA_DHTconv = zeros(1,N);
RM_FFTconv = zeros(1,N);  RA_FFTconv = zeros(1,N);
RM_BFconv = zeros(1,N);   RA_BFconv = zeros(1,N);

for i = 1 : N            % The cases of x   
    K = Lh + J - 1;      % K is the the block length  
    L = ceil(Lx(i)/J);   % Parameter L, L is the number of segements    

% RV-based convolution
    % (optimal)Combine FFT and Matveev's,Recursively break down DFT into smaller DFTs.
    % Use Matveev's eigenvectors for smaller DFTs. RV_FFT
    inp = randi([-1 1]*127, K, 1);    
    myfftpfa(inp(1: K));
    RM_ORRVconv(i,1) = 2*L*a(1)+4*L*K;
    RA_ORRVconv(i,1) = 2*L*a(2)+2*L*K;
    a = [0;0];  

    % Matveev's eigenvectors -> y(:,2)  RV_DFT
    ops = numops(K,'r');
    RM_RVconv(i,1) = 2*L*ops(1)+4*L*K;
    RA_RVconv(i,1) = 2*L*ops(2)+2*L*K;

% Optimal DHT-based convolution    
    % Make K a power of 2
    K = 2^(ceil(log2(K)));

    % Mutiplication number and Addition number
    RM_DHTconv(i,1) =  2*L*(2*K*log2(K)-6*K+8)+ 4*L*K;
    RA_DHTconv(i,1)  = 2*L*(3*K*log2(K)-3*K+4)+ 2*L*K;    
     
% FFT-based convolution
    N(i) = Lx(i) + Lh - 1;
    % Make N to a power of 2
    N(i) = 2^(ceil(log2(N(i))));

    % Radix 2 FFT Mutiplication number and Addition number
    RM_FFTconv(i, 1) = 2*(2*N(i)*log2(N(i))-6*N(i)+8)+4*N(i);  
    RA_FFTconv(i, 1) = 2*(3*N(i)*log2(N(i))-3*N(i)+4)+2*N(i);   

% Brute force convolution (x*convmtx(h,Lx))
    N(i) = Lx(i) + Lh - 1;
    if Lh < Lx(i)
        RM_BFconv(i,1) = Lh^2 + Lh*(Lx(i)-Lh-1);       
        RA_BFconv(i,1) = (Lh-1)^2 + (Lh-1)*(Lx(i)-Lh-1); 
    else
        RM_BFconv(i,1) = 2*sum(1:Lx(i))+ Lx(i)*(Lh-Lx(i)-1); 
        RA_BFconv(i,1) = 2*sum(1:Lx(i)-1)+ (Lx(i)-1)*(Lh-Lx(i)-1);        
    end
    
end

%% Display 
figure('color', [1 1 1]);
%plot(Lx+Lh-1, RM_BFconv(:,1)); hold on;
plot(Lx+Lh-1, RM_FFTconv(:,1)); hold on;
plot(Lx+Lh-1, RM_DHTconv(:,1)); hold on;
plot(Lx+Lh-1, RM_ORRVconv(:,1)); hold on;
plot(Lx+Lh-1, RM_RVconv(:,1)); hold on;
legend('Cooley-Tukey FFT-based Conv','DHT-based Conv','RV-FFT based Conv','RV-DFT based Conv');
set(legend,...
    'Position',[0.28242284406258 0.785834710743801 0.109392492400644 0.108550660519981]);
xlabel('Convolution length'); ylabel('No. of Multiplications');
title('Lh=12 ,J=13,vary Lx');

figure('color', [1 1 1]);
%plot(Lx+Lh-1, RM_BFconv(:,1)); hold on;
plot(Lx+Lh-1, RA_FFTconv(:,1)); hold on;
plot(Lx+Lh-1, RA_DHTconv(:,1)); hold on;
plot(Lx+Lh-1, RA_ORRVconv(:,1)); hold on;
plot(Lx+Lh-1, RA_RVconv(:,1)); hold on;
legend('Cooley-Tukey FFT-based Conv','DHT-based Conv','RV-FFT based Conv','RV-DFT based Conv');
set(legend,...
    'Position',[0.28242284406258 0.785834710743801 0.109392492400644 0.108550660519981]);
xlabel('Convolution length'); ylabel('No. of Additions');
title('Lh=12,J=13,vary Lx');