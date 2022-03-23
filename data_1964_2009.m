% 
% Create detrended time series for technology and unemployment using US data for 1964:Q1–2009:Q2
% ux: detrended unemployment series, with mean 5.8%
% ax: detrended technology series, with mean 1
% amc: AX cast into the state space A for the technology Markov Chain
% 

global nsample rho_a sigma_a

nsample=182;
[rho_a,sigma_a,ax]=TECHNO(0.666,nsample);

%% Productivity

% Interpolate weekly series from quarterly series to simulate model
del=ax(2:end)-ax(1:end-1);
del=del'./12;
DEL=[zeros(size(del));del; 2.*del; 3.*del; 4.*del; 5.*del; 6.*del; 7.*del; 8.*del; 9.*del; 10.*del; 11.*del];
AX=repmat(ax(1:end-1)',12,1);
WAX=AX+DEL;
WAX=WAX(1:end);
wax=[WAX,ax(end)];
wnsample=size(wax,2);

% Stochastic process replicating labor productivity at weekly frequency 
ns=200; % Number of states in Markov chain 
[Z,PI,EPS]=MAKEMC(ns,rho_a,sigma_a);
A=exp(Z);

%% Construct technology in space of Markov chain

amc=[]; % Quarterly series
for i=1:max(size(ax))
	res=find(ax(i)<A',1,'first');
	amc=[amc,res];
end

wamc=[]; % Weekly series
for i=1:max(size(wax))
	res=find(wax(i)<A',1,'first');
	wamc=[wamc,res];
end

%% Unemployment

% Get unemployment data (monthly, from BLS)  
fid=fopen('data/CPS-UR.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
UR=TX{4}; % Unemployment rate
UR=QUARTER(UR); % Make quarterly averages
UR=log(UR);
UR=UR(end-(nsample)+1:end);

uhp=hpfilter(UR,10^5); % Detrended unemployment
ux0=exp(uhp);
NAIRU=exp(UR)./ux0; % Trend
nairu=mean(NAIRU)./100 % Check that it is at 5.8%
ux=ux0.*nairu; % Detrended unemployment -- mean is average unemployment


% Interpolate weekly series from quarterly series 
del=ux(2:end)-ux(1:end-1);
del=del'./12;
DEL=[zeros(size(del));del; 2.*del; 3.*del; 4.*del; 5.*del; 6.*del; 7.*del; 8.*del; 9.*del; 10.*del; 11.*del];
UX=repmat(ux(1:end-1)',12,1);
WUX=UX+DEL;
WUX=WUX(1:end);
wux=[WUX,ux(end)];