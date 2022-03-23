% 
% Equilibrium diagrams for various search-and-matching models
% 

% General setup
clear all;close all;
setup_equilibrium;

UG=[0.01:0.001:0.18];
NG=sort((1-UG)./(1-s));
THG=uinv(UG);

% Model with wage rigidity
KMPS=r.*cMPS./q(THG);
JMPS=(a-wMPS).*ones(size(UG));

% Job-rationing model
KJR=r.*c./q(THG);
JJR=B.*(1-UG).^(alpha-1) -w;

% Canonical model
KMP=cMP./q(THG).*r+cMP.*betaMP.*(1-s).*delta.*THG
JMP=(1-betaMP).*ones(size(UG));

% Model with diminishing returns
JSZ=(1-betaSZ)./(1-betaSZ.*(1-alphaSZ)).*alphaSZ.*(1-UG).^(alphaSZ-1)./(1-s).^(alphaSZ-1);
KSZ=cSZ./q(THG).*r+cSZ.*betaMP.*(1-s).*delta.*THG;

%% Plot equilibrium diagrams

% Model with wage rigidity
figure(2)
clf
plot(NG,JMPS(end:-1:1),'b','LineWidth',4)
hold on
plot(NG,KMPS(end:-1:1),'--r','LineWidth',4)
xlabel('Employment','FontSize',22)
ylabel('Model with wage rigidity','FontSize',22)
ylim([0,.05])
xlim([0.85,1])
set(gca,'YGrid','on')
set(gca,'XGrid','on')
set(gca,'FontSize',22)
h_legend=legend('Gross marginal profit','Marginal recruiting expenses');
set(h_legend,'FontSize',22,'Location','NorthWest');  
print('-depsc','MPSAERq.eps')

% Canonical model
figure(1)
clf
plot(NG,JMP(end:-1:1),'b','LineWidth',4)
hold on
plot(NG,KMP(end:-1:1),'--r','LineWidth',4)
xlabel('Employment','FontSize',22)
ylabel('Canonical model','FontSize',22)
ylim([0,.2])
xlim([0.85,1])
set(gca,'YGrid','on')
set(gca,'XGrid','on')
set(gca,'FontSize',22)
h_legend=legend('Gross marginal profit','Marginal recruiting expenses');
set(h_legend,'FontSize',22,'Location','NorthWest');  
print('-depsc','MPAERq.eps')

% Canonical model with lower recruiting costs
figure(3)
clf
plot(NG,JMP(end:-1:1),'b','LineWidth',4)
hold on
plot(NG,0.15.* KMP(end:-1:1),'--r','LineWidth',4)
plot(NG,0.5.* KMP(end:-1:1),':r','LineWidth',4)
plot(NG,KMP(end:-1:1),':r','LineWidth',4)
xlabel('Employment','FontSize',22)
ylabel('Canonical model','FontSize',22)
ylim([0,.2])
xlim([0.85,1])
set(gca,'YGrid','on')
set(gca,'XGrid','on')
set(gca,'FontSize',22)
h_legend=legend('Gross marginal profit','Marginal recruiting expenses');
set(h_legend,'FontSize',22,'Location','NorthWest');  
print('-depsc','MP2AERq.eps')

% Model with diminishing returns
figure(4)
clf
plot(NG,JSZ(end:-1:1),'b','LineWidth',4)
hold on
plot(NG,KSZ(end:-1:1),'--r','LineWidth',4)
xlabel('Employment','FontSize',22)
ylabel('Model with diminishing returns','FontSize',22)
ylim([0,.2])
xlim([0.85,1])
set(gca,'YGrid','on')
set(gca,'XGrid','on')
set(gca,'FontSize',22)
h_legend=legend('Gross marginal profit','Marginal recruiting expenses');
set(h_legend,'FontSize',22,'Location','NorthWest');  
print('-depsc','SZAERq.eps')

% Job-rationing model
figure(11)
clf
plot(NG,JJR(end:-1:1)+0.003,'b','LineWidth',4)
hold on
plot(NG,KJR(end:-1:1),'--r','LineWidth',4)
xlabel('Employment','FontSize',22)
ylabel('Model with job rationing','FontSize',22)
ylim([0,.05])
xlim([0.85,1])
set(gca,'YGrid','on')
set(gca,'XGrid','on')
set(gca,'FontSize',22)
h_legend=legend('Gross marginal profit','Marginal recruiting expenses');
set(h_legend,'FontSize',22,'Location','NorthWest');  
print('-depsc','JRAERq1.eps')

%% Plot recession in job-rationing model

uz=0.1;
thz=uinv(uz);
res0=1;
OBJF=@(ax)(r.*c.*ax./q(thz)-(ax.*B.*(1-uz).^(alpha-1)-w.*ax.^gamma)).^2;
[res,val,exitflag]=fsolve(OBJF,res0,optimset('TolFun',10^(-16)));
az=res
JJR2=az.*B.*(1-UG).^(alpha-1) -w.*az^(gamma);

figure(12)
clf
plot(NG,JJR2(end:-1:1),'b','LineWidth',4)
hold on
plot(NG,KJR(end:-1:1),'--r','LineWidth',4)
plot(NG,JJR(end:-1:1)+0.003,':b','LineWidth',4)
xlabel('Employment rate','FontSize',22)
ylabel('Model with job rationing','FontSize',22)
ylim([0,.05])
xlim([0.85,1])
set(gca,'YGrid','on')
set(gca,'XGrid','on')
set(gca,'FontSize',22)
h_legend=legend('Gross marginal profit','Marginal recruiting expenses');
set(h_legend,'FontSize',22,'Location','NorthWest');  
print('-depsc','JRAERq2.eps')