% 
% Moments of US quarterly data
% Period: 1964:Q1–2009:Q2
% HP-filtered with weight 1600
% 

clear all;close all;

whp=1600; % HP-filter weight, standard in literature for quarterly data

%% Get nominal wage data (CES, monthly, 1964:M1–2009:M6)

fid=fopen('data/CES-HWAGEPROD.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
NW=TX{4};
NW=QUARTER(NW); % Make quarterly averages

%% Get price data (CPI, monthly, 1947:M1–2009:M6)

fid=fopen('data/CPI-URBAN.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
CPI=TX{4};
CPI=QUARTER(CPI);

%% Construct (log) real wage series

ind=max(size(NW));
CPI=CPI(end-ind+1:end);
RW=log(NW)-log(CPI);

%% Vacancies (merge two datasets: Conference Board and JOLTS)

% Get vacancy data from Conference Board (index, monthly, 1951:M1–2006:M7)  
fid=fopen('data/HELPWANT.txt');
TX = textscan(fid,'%f %f %f %f','HeaderLines',0,'delimiter', '\t');
fclose(fid);
V=TX{4};
V=QUARTER(V); % Make quarterly averages: 1951:Q1–2006:Q2
Vtest=log(V(end-22+1:end));
V=log(V(1:end-22));

% Get JOLTS (monthly, 1,000, 2000:M12–2009:M6)
fid=fopen('data/JOLTS-JOLNF.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
VJ=TX{4}; % Vacancies
VJ=QUARTER(VJ(2:end)); % Make quarterly averages (get rid of 2000:M12)
VJtest=log(VJ(1:22));VJ=log(VJ);

V=[V;VJ-VJtest(1)+Vtest(1)]; % Scale JOLTS series

%% Get CPS data (1948:M1–2009:M6)

fid=fopen('data/CPS-UL.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
UL=TX{4}; % Unemployment level
UL=QUARTER(UL); % Make quarterly averages
UL=log(UL);

fid=fopen('data/CPS-UR.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
UR=TX{4}; % Unemployment rate
UR=QUARTER(UR); % Make quarterly averages
UR=log(UR);

%% Labor market tightness

% Adjust size
ind=max(size(V));
UL=UL(end-ind+1:end);
% Combine series
TH=V-UL; % In logs

%% Productivity data (quarterly, from MSPC BLS, 1947:Q1–2010:Q2)

fid=fopen('data/MSPC-OUTPUT.txt');
TX= textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', '\t');
fclose(fid);
Y=TX{4};

fid=fopen('data/MSPC-EMP.txt');
TX= textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', '\t');
fclose(fid);
N=TX{4};

alpha=0.666; % In the job-rationing model
Y=log(Y);
N=log(N);
A=Y-alpha.*N;

%% Adjust size of series, log, and HP filter

ind=max(size(RW)); % Nominal wage series is shortest
A=A(end-ind+1:end);
Y=Y(end-ind+1:end);
V=V(end-ind+1:end);
UR=UR(end-ind+1:end);
TH=TH(end-ind+1:end);

% Detrended with HP filter
a=hpfilter(A,whp);
rw=hpfilter(RW,whp);
unemp=hpfilter(UR,whp);
y=hpfilter(Y,whp);
th=hpfilter(TH,whp);
v=hpfilter(V,whp);

% Tabulate moments for 1964:Q1–2009:Q2
DD=[unemp,v,th,rw,y,a];
[moy,dev,autoc,Q]=SUMSTAT(DD);
TAB3=[dev',autoc',Q];
fid = fopen('moments_1600.txt', 'wt');
fprintf(fid, '& %4.3f  & %4.3f & %4.3f  & %4.3f & %4.3f & %4.3f \\\\ \n', TAB3);
fclose(fid);