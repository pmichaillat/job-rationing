%==========================================================
% Moments of US quarterly data, HP(10^5)-filtered
% focus on U,V,TH,W,Y,A 
% period: 1964:Q1--2009:Q2
%==========================================================

clear all;close all;

whp=1600; %more standard in BC literature for quarterly data

%%===================================================================   get nominal wage data (CES, monthly)  01/1964--06/2009
fid=fopen('data/CES-HWAGEPROD.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
NW=TX{4};
NW=QUARTER(NW); %make quarterly averages

%%===================================================================   get price data (CPI, monthly)  01/1947--06/2009
fid=fopen('data/CPI-URBAN.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
CPI=TX{4};
CPI=QUARTER(CPI);

%%===============================================     Construct (log) real wage series
ind=max(size(NW));
CPI=CPI(end-ind+1:end);
RW=log(NW)-log(CPI);

%%===============     Vacancies (merge two datasets: conference board and JOLTS)
%%%%%%%% get vacancy data from Conference Board(index, monthly, 01/1951--07/2006)  %%%%%%%%%%%
fid=fopen('data/HELPWANT.txt');
TX = textscan(fid,'%f %f %f %f','HeaderLines',0,'delimiter', '\t');
fclose(fid);
V=TX{4};
V=QUARTER(V); %make quarterly averages: 1951:Q1--2006:Q2
Vtest=log(V(end-22+1:end));
V=log(V(1:end-22));

%%%%%%%%%%%%  get JOLTS (monthly, 1,000, 12/2000 -- 06/2009) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('data/JOLTS-JOLNF.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
VJ=TX{4};%vacancies
VJ=QUARTER(VJ(2:end)); %make quarterly averages -- get rid of 12/2000
VJtest=log(VJ(1:22));VJ=log(VJ);

V=[V;VJ-VJtest(1)+Vtest(1)];%scale jolts series


%%===================================================================   get CPS data 01/1948 -- 06/2009
fid=fopen('data/CPS-UL.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
UL=TX{4}; %unemployment level
UL=QUARTER(UL);%make quarterly averages
UL=log(UL);

fid=fopen('data/CPS-UR.txt');
TX = textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', ',');
fclose(fid);
UR=TX{4}; %unemployment rate
UR=QUARTER(UR);%make quarterly averages
UR=log(UR);

%%===================================================================    Labor market tightness
%adjust size
ind=max(size(V));
UL=UL(end-ind+1:end);
%combine series
TH=V-UL;%in logs - ratio of levels


%%===================================================================     productivity data (quarterly, from MSPC BLS, 1947:Q1--2010:Q2)
fid=fopen('data/MSPC-OUTPUT.txt');
TX= textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', '\t');
fclose(fid);
Y=TX{4};

fid=fopen('data/MSPC-EMP.txt');
TX= textscan(fid,'%s %s %s %f','HeaderLines',1,'delimiter', '\t');
fclose(fid);
N=TX{4};

alpha=0.666;%in the Michaillat(2010) model
Y=log(Y);
N=log(N);
A=Y-alpha.*N;%could also compute Y/N if required

%%======================================     Adjust size of series, log, and HP filter
ind=max(size(RW));%nominal wage series is shortest
A=A(end-ind+1:end);
Y=Y(end-ind+1:end);
V=V(end-ind+1:end);
UR=UR(end-ind+1:end);
TH=TH(end-ind+1:end);

%detrended with HP filter
a=hpfilter(A,whp);
rw=hpfilter(RW,whp);
unemp=hpfilter(UR,whp);
y=hpfilter(Y,whp);
th=hpfilter(TH,whp);
v=hpfilter(V,whp);

%tabulate moments for 1964:Q1--2009:Q2
DD=[unemp,v,th,rw,y,a];
[moy,dev,autoc,Q]=SUMSTAT(DD);
TAB3=[dev',autoc',Q];
fid = fopen('table/US_moments_1600.txt', 'wt');
fprintf(fid, '& %4.3f  & %4.3f & %4.3f  & %4.3f & %4.3f & %4.3f \\\\ \n', TAB3);
fclose(fid);


