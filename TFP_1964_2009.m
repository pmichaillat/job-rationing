%%==============================================================================
%% Create detrended time series for TFP and unemployment using US data for 1964:Q1--2009:Q2
%% ux: detrended unemployment series, with mean 5.8%
%% ax: detrended TFP series, with mean 1
%% amc: AX cast into the state space A for the technology Markov Chain
%%==============================================================================

global nsample

nsample=182;

%%=========================================================================   Productivity

%%%%%%%%%%%  Get TFP data (quarterly, from John Fernald)     %%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen('/Users/Pascal/Documents/Data/US_DATA/FERNALD.txt');
BB= textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',1,'delimiter', '\t');
fclose(fid);
BB=BB(2:end);%supress dates
TFP=BB{9};%400*dln(TFP)
TFP=TFP./400;%
TFP=cumsum(TFP);%log(TFP)
AA=TFP(end-(182)+1:end);%period of interest (goes to 2009:Q2)
ahp=hpfilter(AA,10^5);%detrended TFP
ax=exp(ahp);%mean of 1, as normalized


%%%%%%%%%%%  Interpolate weekly series from quarterly series to simulate model   %%%%%%%%%%%%%%%%%%%%%%%%%%%

del=ax(2:end)-ax(1:end-1);
del=del'./12;
DEL=[zeros(size(del));del; 2.*del; 3.*del; 4.*del; 5.*del; 6.*del; 7.*del; 8.*del; 9.*del ; 10.*del; 11.*del];%one week= 1/12 quarter
AX=repmat(ax(1:end-1)',12,1);
WAX=AX+DEL;
WAX=WAX(1:end);
wax=[WAX,ax(end)];

wnsample=size(wax,2);

%%%%%%%%%%%  Run AR(1) regression on detrended technology series   %%%%%%%%%%%%%%%%%%%%%%%%%%%
[rho,sigma]=AR(ahp,1)

%%%%%%%%%%%   stochastic process -- replicates AL at weekly frequency  %%%%%%%%%%%%%%%%%%%%%%%%%%

ns=200;%number of states in MC 
rho_a=(rho).^(1./12) ;  % persistence in labor productivity (estimated for 1964-2009 with AR(1) and brought at weekly)
sigma_a=sigma./(1+rho_a^2+rho_a^4+rho_a^6+rho_a^8+rho_a^10+rho_a^12+rho_a^14+rho_a^16+rho_a^18+rho_a^20+rho_a^22).^(0.5) ; %variance of epsilon in log labor productivity (HM08b)
sqrt(sigma_a^2./(1-rho_a^2));
[Z,PI,EPS]=MAKEMC(ns,rho_a,sigma_a);
A=exp(Z);

%%===================================================================    Unemployment

%%%%%%%%%%%  Get unemployment data  (monthly, from BLS)     %%%%%%%%%%%%%%%%%%%%%%%%%%%


fid=fopen('/Users/Pascal/Documents/Data/US_DATA/PC.txt');
BB = textscan(fid,'%s %s %f %f %f','HeaderLines',0,'delimiter', '\t');
fclose(fid);
BB=QUARTER(BB(5),0);
UNF=BB{1}./100;
UNF=UNF(end-(nsample)+1:end);
pas=max(size(UNF));
ulog=log(UNF);
uhp=hpfilter(ulog,10^5);%detrended unemp.
ux0=exp(uhp);
NAIRU=UNF./ux0;%trend
nairu=mean(NAIRU) %check that it is at 5.8%
ux=ux0.*nairu;%detrended unemployment -- mean is average unemployment (nairu)

%%%%%%%%%%%  Interpolate weekly series from quarterly series  %%%%%%%%%%%%%%%%%%%%%%%%%%%

del=ux(2:end)-ux(1:end-1);
del=del'./12;
DEL=[zeros(size(del));del; 2.*del; 3.*del; 4.*del; 5.*del; 6.*del; 7.*del; 8.*del; 9.*del ; 10.*del; 11.*del];%one week= 1/12 quarter
UX=repmat(ux(1:end-1)',12,1);
WUX=UX+DEL;
WUX=WUX(1:end);
wux=[WUX,ux(end)];


%%==============================================     Construct A  in space of MC

amc=[];%quarterly series
for i=1:max(size(ax))
	res=find(ax(i)<A',1,'first');
	amc=[amc,res];
end

wamc=[];%weekly series
for i=1:max(size(wax))
	res=find(wax(i)<A',1,'first');
	wamc=[wamc,res];
end