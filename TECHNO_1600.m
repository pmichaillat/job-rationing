% 
% Construct the time series for technology from US data, HP-filter it, and report detrended series as well as moments of series 
% Moments are obtained by estimating a weekly AR(1) process
% Focus on the period 1964:Q1–2009:Q2
% HP-filter weight is 1600
% 

function [rho_a,sigma_a,ax]=TECHNO_1600(alpha,nsample)

whp=1600; % Weight on HP filter

%% Get productivity data (quarterly, from BLS)

fid=fopen('data/MSPC-OUTPUT.txt');
TX= textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', '\t');
fclose(fid);
Y=TX{4};

fid=fopen('data/MSPC-EMP.txt');
TX= textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', '\t');
fclose(fid);
N=TX{4};

Y=log(Y);
N=log(N);
A=Y-alpha.*N;

A=A(end-(nsample)+1:end);
ahp=hpfilter(A,whp); % Detrended log technology
ax=exp(ahp); % Detrended technology: mean of 1, as normalized

%% Run AR(1) regression on detrended technology series 

[rho,sigma]=AR(ahp,1);
rho_a=(rho).^(1./12); % Persistence in labor productivity (estimated for 1964–2009 with AR(1) and brought at weekly frequency)
sigma_a=sigma./(1+rho_a^2+rho_a^4+rho_a^6+rho_a^8+rho_a^10+rho_a^12+rho_a^14+rho_a^16+rho_a^18+rho_a^20+rho_a^22).^(0.5); % Variance of epsilon in log labor productivity