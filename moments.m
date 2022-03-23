% 
% moments of linear DSGE model (already in logs)
% 

clear all;close all;
setupsimul;
whp=10^5; % Weight on HP filter from Shimer (2005)
rep=30; % Number of replications
cut=100;
T=rep.*182*(3*4)+cut*(3*4) % Samples of 182 quarters
Eps=sigma_a.*randn(1,T); % Realization of errors

%% Simulation of weekly time series

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

%% Make quarterly values

tht=tht(1:4:end);
ut=ut(1:4:end);
vt=vt(1:4:end);

% Quarterly averages
tht=1./3.*(tht(1:3:end-2)+tht(2:3:end-1)+tht(3:3:end));
ut=1./3.*(ut(1:3:end-2)+ut(2:3:end-1)+ut(3:3:end));
vt=1./3.*(vt(1:3:end-2)+vt(2:3:end-1)+vt(3:3:end));

at=at(1:12:end);
wt=wt(1:12:end);
yt=yt(1:12:end);

%% Moments as averages over samples of 182 quarters

ix=0;moy=[];dev=[];autoc=[];Q=[];

for i=1:rep
ran=[1+ix:182+ix];
D=[ut(ran)',vt(ran)',tht(ran)',wt(ran)',yt(ran)',at(ran)'];
% HP-filter all samples of model-generated series
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
fid = fopen('mean.txt', 'wt');
fprintf(fid, '& %4.3f  & %4.3f & %4.3f & %4.3f  & %4.3f & %4.3f \\\\ \n', TAB3);
fclose(fid);

TAB3=[dev2,autoc2,Q2];
fid = fopen('variance.txt', 'wt');
fprintf(fid, ' & (%4.3f) & (%4.3f) & (%4.3f) & (%4.3f)  & (%4.3f) & (%4.3f)\\\\ \n', TAB3);
fclose(fid);