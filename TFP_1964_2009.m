% 
% Create detrended time series for TFP and unemployment using US data for 1964:Q1–2009:Q2
% ux: detrended unemployment series, with mean 5.8%
% ax: detrended TFP series, with mean 1
% amc: AX cast into the state space A for the technology Markov Chain
% Use utilization-adjusted TFP data constructed by Fernald (2009)
% 

global nsample

nsample=182;

%% Technology

% Get utilization-adjusted TFP data from Fernald (2009)
fid=fopen('data/FERNALD.txt');
BB= textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',1,'delimiter', '\t');
fclose(fid);
BB=BB(2:end); % Suppress dates
TFP=BB{9}; % 400*dln(TFP)
TFP=TFP./400;
TFP=cumsum(TFP); % log(TFP)
AA=TFP(end-(182)+1:end); % Period of interest (goes to 2009:Q2)
ahp=hpfilter(AA,10^5); % Detrended TFP
ax=exp(ahp); % Mean of 1, as normalized

% Interpolate weekly series from quarterly series to simulate model   
del=ax(2:end)-ax(1:end-1);
del=del'./12;
DEL=[zeros(size(del));del; 2.*del; 3.*del; 4.*del; 5.*del; 6.*del; 7.*del; 8.*del; 9.*del; 10.*del; 11.*del];
AX=repmat(ax(1:end-1)',12,1);
WAX=AX+DEL;
WAX=WAX(1:end);
wax=[WAX,ax(end)];
wnsample=size(wax,2);

% Run AR(1) regression on detrended technology series for 1964–2009 
[rho,sigma]=AR(ahp,1);

% Stochastic process that replicates labor productivity at weekly frequency
ns=200; % Number of states in Markov chain 
rho_a=(rho).^(1./12); % Persistence in labor productivity, weekly frequency
sigma_a=sigma./(1+rho_a^2+rho_a^4+rho_a^6+rho_a^8+rho_a^10+rho_a^12+rho_a^14+rho_a^16+rho_a^18+rho_a^20+rho_a^22).^(0.5); % Variance of epsilon in log labor productivity
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
UNF=TX{4}./100; % Unemployment rate
UNF=QUARTER(UNF); % Make quarterly averages
UNF=UNF(end-(nsample)+1:end);
pas=max(size(UNF));
ulog=log(UNF);
uhp=hpfilter(ulog,10^5); % Detrended unemployment
ux0=exp(uhp);
NAIRU=UNF./ux0; % Trend
nairu=mean(NAIRU);
ux=ux0.*nairu; % Detrended unemployment -- mean is average unemployment

% Interpolate weekly series from quarterly series  
del=ux(2:end)-ux(1:end-1);
del=del'./12;
DEL=[zeros(size(del));del; 2.*del; 3.*del; 4.*del; 5.*del; 6.*del; 7.*del; 8.*del; 9.*del; 10.*del; 11.*del];
UX=repmat(ux(1:end-1)',12,1);
WUX=UX+DEL;
WUX=WUX(1:end);
wux=[WUX,ux(end)];