%% FFTs_num_mults_Fig1

% 7/17/2023 FAMU-FSU colledge of engineering
% This file is programmed by Weiwei Wang (ww20br@fsu.edu).
% Supervised by Prof. Victor DeBrunner (victor.debrunner@eng.famu.fsu.edu).
% Reviewed by Prof. Linda DeBrunner (linda.debrunner@eng.famu.fsu.edu).

% This code is used to do the number of multiplications comparison for 
% different DFTs.

%%
clear
close all

lb = 3; ub = 64;  % lower and upper bounds
n = (lb:ub)';
%m = floor([n+4, n+2, n-1, n+1]/4);

y = zeros(length(n),4);

% Combine Prime Factor Algorithm (PFA) and Matveev's -> y(:,1)
% Recursively break down DFT into smaller DFTs.
% Use Matveev's eigenvectors for smaller DFTs.
clear global a;
global a;
inp = randi([-1 1]*127, ub, 1);
for i = n'
    myfftpfa(inp(1:i));
    y(i-lb+1,1) = a(1);
    a = [0;0];
end

% Matveev's eigenvectors -> y(:,2)
for i = n'
    ops = numops(i,'r');
    y(i-lb+1,2) = ops(1);
end

% Brute force -> y(:,3)
y(:,3) = 2*(n-1).^2;

% Cooley-Tukey FFT -> y(:,4)
% for length N, no. of mult = 2N*log2(N) - 6N + 8
lb2 = ceil(log2(lb));
ub2 = floor(log2(ub));
n2 = 2.^(lb2:ub2)';
y(:,4) = NaN;
for i = n2'
    y(i-lb+1,4) = 2*i*log2(i) - 6*i + 8;
end

% f = @(xaxis, yaxis) createfigure(xaxis, yaxis);
createfigure(n, y,lb,ub);
% plot(n,y,'ro','bo','mx','kx','g*');
clearvars lb ub m lb2 ub2 n2 inp i


%%  createfigure(X1, YMatrix1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

function createfigure(X1, YMatrix1,lb,ub)

% Create figure
figure1 = figure('OuterPosition',[720.2 301 516.8 456.8]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineStyle','none');
set(plot1(1),'Marker','*','Color',[0.266 0.95 0.188]);% Combine FFT and Matveev's -> y(:,1)% Recursively break down DFT into smaller DFTs.
set(plot1(2),'Marker','o','Color',[0.9 0.325 0.098]);% Matveev's eigenvectors -> y(:,2)
set(plot1(3),'Marker','x','Color',[0.8 0.6 0.1]);% Brute force -> y(:,5)
set(plot1(4),'Marker','diamond','Color',[0 0 0]);% Cooley-Tukey FFT -> y(:,6)

% set legend
leg = {'RV-FFT', 'RV-DFT', 'Brute force DFT','Cooley-Tukey FFT'};
legend([plot1(3) plot1(2) plot1(1) plot1(4)],{leg{3},leg{2}, leg{1}, leg{4}})

% Create ylabel
ylabel('No. of multiplications');

% Create xlabel
xlabel('Size of DFT');

% Create title
title('Number of multiplications');

% Set axes limits
axis([X1(1) X1(end) 0 max(YMatrix1(:)+100)]);

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.14242284406258 0.679834710743801 0.309392492400644 0.208550660519981]);
end
