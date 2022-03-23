% 
% Simulate with shooting the path of unemployment and labor market tightness
% Shocks to technology are actual shocks from US data -- Markov Chain (transition matrix PI) 
% TH: solution of stochastic equilibrium
% A: discrete technology process
% PI: Markov transition matrix
% TH: stochastic equilibrium (used for 1st guess)
% T: length of simulation
% amc: time series of states
% n0: employment in past period
% output: time series for unemployment ut, labor market tightness tht, technology at, measure of layoffs sigmat
% 

function [Yt]=SIMULFT(TH,A,wamc,n0,EY,YLR)

global delta eta c a alpha sigma_a omega rho_a s markup r varsigma  u_target th_target sigma w gamma
global q f u finv qinv uinv
global ns ynum
global apos thpos npos mplpos hpos wpos Rpos  upos


%% Key parameters

k0=25; % Number of expectations to compute to have stable result
hori=1;
k1=k0+2.*hori+1; % Total number of expectations to compute (including current period)

n=max(size(wamc))

%% Initialization

Yt=zeros(ynum,n+1);
Y0=zeros(ynum,1); % Vector Y at t-1 (only the value of employment matter)
Y0(npos)=n0;
Yt(:,1)=Y0;
s1=wamc(1); % First state
kIII=k1+1;

%% Iteration over 1964–2009 period to solve model

for i=1:n 
	si=wamc(i);
	grt=squeeze(EY(:,si,:)); % Include guess for current period
  	[Yi,er,kIII,eII,eIII]=SHOOTING(grt,Yt(:,i),si,floor(0.5.*k1+0.5.*(kIII-1))); % Include type I, II, and III iteration
	Yt(:,i+1)=Yi;
	figure(4)
	scatter([i],[Yi(upos)],'red')
	eII
	eIII
	hold on
end

Yt=Yt(:,2:end); % Get rid of initial value Y0