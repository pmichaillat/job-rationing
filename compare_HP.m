%%==============================================================================
%% Redo analysis with HP weight of 1600
%% Checks that moments are similar, fit of model (comparison of moments), simulated and actual unemployment series
%%==============================================================================

clear all;close all;

%%==============================================================================
%%                                                                  setupsimul.m
%%==============================================================================

global  W_bar C_bar TH_bar  N_bar H_bar U_bar UC_bar UF_bar Y_bar 
global rho_a sigma_a nsample
global wpos cpos thpos npos hpos upos ypos apos  a_pos

%%==============================================================     Calibration

setup;%model parameters
[W_bar,C_bar,TH_bar,N_bar,H_bar,U_bar,UC_bar,UF_bar,Y_bar]=STEADYLL(w,gamma);%steady states
nsample=182;


whp=1600;%weight on hp filter (quarterly frequencies)- shimer (2005): 10^5 - conventional: 1600

%%==============================================================================
%%                                                                      TECHNO.m
%%==============================================================================

%%%%%%%%%%%  Get productivity data (quarterly, from BLS)     %%%%%%%%%%%%%%%%%%%%%%%%%%%

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

A=A(end-(nsample)+1:end);
ahp=hpfilter(A,whp);%detrended log techno
ax=exp(ahp);%detrended techno -- mean of 1, as normalized

%%%%%%%%%%%  Run AR(1) regression on detrended technology series   %%%%%%%%%%%%%%%%%%%%%%%%%%%
[rho,sigma]=AR(ahp,1);
rho_a=(rho).^(1./12)  % persistence in labor productivity (estimated for 1964-2009 with AR(1) and brought at weekly)
sigma_a=sigma./(1+rho_a^2+rho_a^4+rho_a^6+rho_a^8+rho_a^10+rho_a^12+rho_a^14+rho_a^16+rho_a^18+rho_a^20+rho_a^22).^(0.5) %variance of epsilon in log labor productivity 



%%==============================================================================
%%                                                                  end TECHNO.m
%%==============================================================================



%positions of variables
wpos	=2;	% wage
cpos   =  3;               % consumption
thpos   =  1;               % market tightness
npos   =  4;               % employment
hpos   =  5;               % hiring. rate
upos   =  6;               % unemp. rate
ypos  =  7;               % output
apos   =  8;               % techno 
a_pos   = 9;               % techno shock

%%===================================================     Solve log-linear model
[AMAT,BMAT,SS,xeq,xvari]=LIN_DSGE(w,gamma,alpha);% to compute and plot impulse response use reduced form solution: x_t=amat*x_{t-1}+b*shock

%%==============================================================================
%%                                                                  end setupsimul.m
%%==============================================================================

rep=100;%number samples - 10,000 in shimer(2005)
cut=100;
T=rep.*182*(3*4)+cut*(3*4)%samples of 182 quarters
Eps=sigma_a.*randn(1,T);%realization of errors

%%===================================================================    simulation of weekly time series
nvar=max(size(AMAT));
x_past=zeros(nvar,1); %start from steady state
shock=zeros(nvar,1);
X=[];
for i=1:T
  shock(a_pos)=Eps(i);
  x_past=AMAT*x_past+BMAT*shock;
  X=[X,x_past];
end

tht=X(thpos,cut*(3*4)+1:end);
ut=X(upos,cut*(3*4)+1:end);
vt=tht+ut;
yt=X(ypos,cut*(3*4)+1:end);
at=X(apos,cut*(3*4)+1:end);
wt=X(wpos,cut*(3*4)+1:end);

%%===================================================================    make quarterly values
tht=tht(1:4:end);
ut=ut(1:4:end);
vt=vt(1:4:end);

%quarterly averages
tht=1./3.*(tht(1:3:end-2)+tht(2:3:end-1)+tht(3:3:end));
ut=1./3.*(ut(1:3:end-2)+ut(2:3:end-1)+ut(3:3:end));
vt=1./3.*(vt(1:3:end-2)+vt(2:3:end-1)+vt(3:3:end));

at=at(1:12:end);
wt=wt(1:12:end);
yt=yt(1:12:end);




%%===================================================================     moments as averages of sample - sample periods of 182 quarters
ix=0;moy=[];dev=[];autoc=[];Q=[];

for i=1:rep
ran=[1+ix:182+ix];
D=[ut(ran)',vt(ran)',tht(ran)',wt(ran)',yt(ran)',at(ran)'];
%%HP-filter all samples of model-generated series
for j=1:6
	D(:,j)=hpfilter(D(:,j),whp);
end
[a1,a2,a3,a4]=SUMSTAT(D);
moy(:,i)=a1;dev(:,i)=a2;autoc(:,i)=a3;Q(:,:,i)=a4;
ix=ix+182;
end

dev1=mean(dev,2);dev2=std(dev,0,2);
autoc1=mean(autoc,2);autoc2=std(autoc,0,2);
Q1=mean(Q,3);Q2=std(Q,0,3);

TAB3=[dev1,autoc1,Q1];
fid = fopen('table/meanHP1600.txt', 'wt');
fprintf(fid, '& %4.3f  & %4.3f & %4.3f & %4.3f  & %4.3f & %4.3f \\\\ \n', TAB3);
fclose(fid);


TAB3=[dev2,autoc2,Q2];
fid = fopen('table/varianceHP1600.txt', 'wt');
fprintf(fid, ' & (%4.3f) & (%4.3f) & (%4.3f) & (%4.3f)  & (%4.3f) & (%4.3f)\\\\ \n', TAB3);
fclose(fid);

%%==============================================================================
%%                                                Add moments for US data with new HP weight
%%==============================================================================


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
UR=TX{4}; %unemployment level
UR=QUARTER(UR);%make quarterly averages
UR=log(UR);

%%===================================================================    Labor market tightness
%adjust size
ind=max(size(V));
UL=UL(end-ind+1:end)
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
fid = fopen('table/US_moments1600.txt', 'wt');
fprintf(fid, '& %4.3f  & %4.3f & %4.3f  & %4.3f & %4.3f & %4.3f \\\\ \n', TAB3);
fclose(fid);

