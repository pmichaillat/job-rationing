% 
% Create detrended time series for technology, unemployment, and productivity using US data for 1964:Q1–2009:Q2
% ux: detrended unemployment series, with mean 5.8%
% ax: detrended technology series, with mean 1
% prodx: detrended technology series, with mean 1
% 

setup;
nsample=182;
whp=10^5;

%% Get unemployment data (monthly, from BLS)

fid=fopen('data/CPS-UR.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
UR=TX{4}; % Unemployment rate
UR=QUARTER(UR); % Make quarterly averages
UR=log(UR);

UR=UR(end-(nsample)+1:end);
uhp=hpfilter(UR,whp); % Detrended unemployment
ux0=exp(uhp);
ux=ux0.*0.058; % Detrended unemployment -- mean is average unemployment

%% Get productivity data (quarterly, from MSPC BLS, 1947:Q1–2009:Q2) 

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
PROD=Y-N;

A=A(end-(nsample)+1:end);
ahp=hpfilter(A,whp); % Detrended technology
ax=exp(ahp); % Mean of 1

PROD=PROD(end-(nsample)+1:end);
prodhp=hpfilter(PROD,whp); % Detrended productivity
prodx0=exp(prodhp);
prodx=prodx0.*(n_target).^(alpha-1); % Detrended techno

Y=Y(end-(nsample)+1:end);
yhp=hpfilter(Y,whp); % Detrended output
yx0=exp(yhp);
yx=yx0.*(n_target).^(alpha); % Detrended output


%% Get real wage data

fid=fopen('data/CES-HWAGEPROD.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
NW=TX{4};
NW=QUARTER(NW); % Make quarterly averages

fid=fopen('data/CPI-URBAN.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
CPI=TX{4};
CPI=QUARTER(CPI);

ind=max(size(NW));
CPI=CPI(end-ind+1:end);
RW=log(NW)-log(CPI);

rwhp=hpfilter(RW,whp); % Detrended real wage
rwx=w.*exp(rwhp); % Mean of wage

fid=fopen('data/ECI.txt'); % Already quarterly, in logs
TX = textscan(fid,'%s %f %f %f %f %f %f %f %f %f','HeaderLines',1,'delimiter', '\t');
fclose(fid);
ECI=TX{5};
ind=find(ECI>0,1);
ECI=ECI(ind:end);
size_eci=size(ECI,1);
CPI=CPI(end-size_eci+1:end);
ECIRW=ECI-log(CPI);

ecirwhp=hpfilter(ECIRW,whp); % Detrended real wage
ecirwx=w.*exp(ecirwhp); % Mean of wage

%% Get labor market tightness data

fid=fopen('data/HELPWANT.txt');
TX = textscan(fid,'%f %f %f %f','HeaderLines',0,'delimiter', '\t');
fclose(fid);
V=TX{4};
V=QUARTER(V); % Make quarterly averages: 1951:Q1--2006:Q2
Vtest=log(V(end-22+1:end));
V=log(V(1:end-22));

fid=fopen('data/JOLTS-JOLNF.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
VJ=TX{4}; % Vacancies
VJ=QUARTER(VJ(2:end)); % Make quarterly averages (get rid of 2000:M12)
VJtest=log(VJ(1:22));VJ=log(VJ);

V=[V;VJ-VJtest(1)+Vtest(1)]; % Scale JOLTS series

fid=fopen('data/CPS-UL.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
UL=TX{4}; % Unemployment level
UL=QUARTER(UL); % Make quarterly averages
UL=log(UL);

% Adjust size
ind=max(size(V));
UL=UL(end-ind+1:end);
% Combine series
TH=V-UL; % In logs - ratio of levels
TH=TH(end-(nsample)+1:end);

thhp=hpfilter(TH,whp); % Detrended tightness
thx=th_target.*exp(thhp); % Mean of th_target

%% Weekly series

% Technology  
anl=QTOW(ax);

% Productivity
prodnl=QTOW(prodx);

% Unemployment 
unl=QTOW(ux);

% Wages
wnl=QTOW(rwx);
eciwnl=QTOW(ecirwx);

% Tightness
thnl=QTOW(thx);

%% Compute variables

u=log(unl);a=log(anl);prod=log(prodnl);
w=log(wnl);eciw=log(eciwnl);
th=log(thnl);nnl=(1-unl(2:end))./(1-s);nnl=[nnl,nnl(end)];n=log(nnl);nx=nnl(1:12:end);