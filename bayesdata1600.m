%%==============================================================================
%% Create detrended time series for technology, unemployment, and productivity using US data for 1964:Q1--2009:Q2
%% ux: detrended unemployment series, with mean 5.8%
%% ax: detrended technology series, with mean 1
%% prodx: detrended technology series, with mean 1
%%==============================================================================
%close all;clear all;
setup
nsample=182;
whp=1600;
w0=w;
%%=================================================================     Get data
%%%%%%%%%%%  Get unemployment data  (monthly, from BLS)     %%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen('data/CPS-UR.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
UR=TX{4}; %unemployment rate
UR=QUARTER(UR);%make quarterly averages
UR=log(UR);


UR=UR(end-(nsample)+1:end);
uhp=hpfilter(UR,whp);%detrended unemp.
ux0=exp(uhp);
ux=ux0.*0.058;%detrended unemployment -- mean is average unemployment (nairu)

%%%%%%%%%%%  Get productivity data (quarterly, from MSPC BLS, 1947:Q1--2009:Q2)   %%%%%%%%%%%%%%%%%%%%%%%%%%%    
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
A=Y-alpha.*N;%could also compute Y/N if required
PROD=Y-N;


A=A(end-(nsample)+1:end);
ahp=hpfilter(A,whp);%detrended unemp.
ax=exp(ahp);% mean of 1

PROD=PROD(end-(nsample)+1:end);
prodhp=hpfilter(PROD,whp);%detrended unemp.
prodx0=exp(prodhp);
prodx=prodx0.*(n_target).^(alpha-1);%detrended techno

Y=Y(end-(nsample)+1:end);
yhp=hpfilter(Y,whp);%detrended unemp.
yx0=exp(yhp);
yx=yx0.*(n_target).^(alpha);%detrended output


%%%%%%%%%%%  Get real wage data   %%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen('data/CES-HWAGEPROD.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
NW=TX{4};
NW=QUARTER(NW); %make quarterly averages

fid=fopen('data/CPI-URBAN.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
CPI=TX{4};
CPI=QUARTER(CPI);

ind=max(size(NW));
CPI=CPI(end-ind+1:end);
RW=log(NW)-log(CPI);

rwhp=hpfilter(RW,whp);%detrended unemp.
rwx=w0.*exp(rwhp);% mean of w0

%%%%%%%%%%%  Get labor market tightness data   %%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen('data/HELPWANT.txt');
TX = textscan(fid,'%f %f %f %f','HeaderLines',0,'delimiter', '\t');
fclose(fid);
V=TX{4};
V=QUARTER(V); %make quarterly averages: 1951:Q1--2006:Q2
Vtest=log(V(end-22+1:end));
V=log(V(1:end-22));

fid=fopen('data/JOLTS-JOLNF.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
VJ=TX{4};%vacancies
VJ=QUARTER(VJ(2:end)); %make quarterly averages -- get rid of 12/2000
VJtest=log(VJ(1:22));VJ=log(VJ);

V=[V;VJ-VJtest(1)+Vtest(1)];%scale jolts series

fid=fopen('data/CPS-UL.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
UL=TX{4}; %unemployment level
UL=QUARTER(UL);%make quarterly averages
UL=log(UL);

%adjust size
ind=max(size(V));
UL=UL(end-ind+1:end);
%combine series
TH=V-UL;%in logs - ratio of levels
TH=TH(end-(nsample)+1:end);

thhp=hpfilter(TH,whp);%detrended unemp.
thx=th_target.*exp(thhp);% mean of w0


%%=========================================================================   Weekly series

%%%%%%%%%%%  Technology   %%%%%%%%%%%%%%%%%%%%%%%%%%%

anl=QTOW(ax);

%%%%%%%%%%%  Productivity %%%%%%%%%%%%%%%%%%%%%%%%%%%
prodnl=QTOW(prodx);

%%%%%%%%%%%  Unemployment  %%%%%%%%%%%%%%%%%%%%%%%%%%%

unl=QTOW(ux);
wnl=QTOW(rwx);
thnl=QTOW(thx);

%%==============================================================     save matrix

u=log(unl);a=log(anl);prod=log(prodnl); w=log(wnl);th=log(thnl);
nnl=(1-unl(2:end))./(1-s);nnl=[nnl,nnl(end)];n=log(nnl);nx=nnl(1:12:end);

% save('logbayesAER1600.mat','prod','n','w','th')