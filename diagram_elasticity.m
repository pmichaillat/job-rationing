% 
% Diagram for elasticities in various search-and-matching models
% 

% General setup
clear all;close all;
setup_elasticity;
UG=[0.01:0.001:0.12];
THG=uinv(UG);
NG=(1-UG)./(1-s);

% Elasticity in job-rationing model
T= alpha.*NG.^(alpha-1).*q(THG)./(r.*c);
EJR=1./(eta+(1-eta).*(1-alpha).*UG.*T);

% Elasticity in model with wage rigidity
EMPR=(1./eta).*ones(size(UG));;

%% Plot elasticities

figure(2)
clf
plot(UG,EJR,'b','LineWidth',4)
hold on
plot(UG,EMPR,'--r','LineWidth',4)
xlabel('Unemployment rate','FontSize',22)
ylabel('|\epsilon^{\theta}_{c}|','FontSize',22)
ylim([0,2.5])
xlim([0.01,0.12])
set(gca,'YGrid','on')
set(gca,'XGrid','on')
set(gca,'FontSize',22)
h_legend=legend('\alpha<1','\alpha=1');
set(h_legend,'FontSize',22);  
print('-depsc','HELAST.eps')