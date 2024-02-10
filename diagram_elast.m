%==================================================================================================
% Diagram for elasticities in equilibrium
%==================================================================================================

%  %general setup
clear all;close all;
setup_elasticity;

UG=[0.01:0.001:0.12];
THG=uinv(UG);
NG=(1-UG)./(1-s);


%JR
AG=(1./w.*(alpha.*NG.^(alpha-1)-r.*c./q(THG))).^(1./(gamma-1));
T= alpha.*NG.^(alpha-1).*q(THG)./(r.*c);
EJR=1./(eta+(1-eta).*(1-alpha).*UG.*T);
ur=@(ax) max(1-((alpha./w).^(1./(1-alpha)).*ax.^((1-gamma)./(1-alpha))),0);
UR=ur(AG);
UF=UG-UR;

%MPS
EMPR=(1./eta).*ones(size(UG));;

%MP
uMP=u_target;thMP=uinv(u_target);nMP=(1-u_target)./(1-s);
Q=betaMP.*delta.*(1-s).*cMP./(1-betaMP);
EMP=1./(eta+(1-eta).*Q.*thMP);
R=betaSZ.*delta.*(1-s).*(1-betaSZ.*(1-alphaSZ)).*cSZ./((1-betaSZ).*alphaSZ);
ESZ=1./(eta+(1-eta).*R.*thMP.*(nMP).^(1-alphaSZ)+(1-alphaSZ).*(1-eta).*uMP);

%-------------------------------------------------------------------------------- 
%GRAPH with elasticities
%figure(2)
%clf
%plot(AG,EJR,'b','LineWidth',4)
%hold on
%plot(AG,EMP,'-.','Color',[0 0.5 0],'LineWidth',4)
%plot(AG,ESZ,'-.m','LineWidth',4)
%plot(AG,EMPR,'--r','LineWidth',4)
%xlabel('Technology','FontSize',22)
%ylabel('|\epsilon^{\theta}_{c}|','FontSize',22)
%ylim([0,3])
%xlim([0.9,1.1])
%set(gca,'YGrid','on')
%set(gca,'XGrid','on')
%set(gca,'FontSize',22)
%h_legend=legend('JR','MP','SZ','MPR')
%set(h_legend,'FontSize',22);  
%print('-depsc','graph/HELAS.eps')
 % 
 % 
 % figure(2)
 % clf
 % plot(UG,EJR,'b','LineWidth',4)
 % hold on
 % plot(UG,EMPR,'--r','LineWidth',4)
 % xlabel('Unemployment rate','FontSize',22)
 % ylabel('|\epsilon^{\theta}_{c}|','FontSize',22)
 % ylim([0,2.5])
 % xlim([0.01,0.12])
 % set(gca,'YGrid','on')
 % set(gca,'XGrid','on')
 % set(gca,'FontSize',22)
 % h_legend=legend('\alpha<1','\alpha=1')
 % set(h_legend,'FontSize',22);  
 % print('-depsc','graph/HELAST.eps')

 figure(3)
 clf
 plot(UG,UF,'b','LineWidth',4)
 hold on
 plot(UG,UG,'--r','LineWidth',4)
 xlabel('Unemployment rate','FontSize',22)
 ylabel('Frictional unemployment u^F','FontSize',22)
 ylim([0.01,0.12])
 xlim([0.01,0.12])
 set(gca,'YGrid','on')
 set(gca,'XGrid','on')
 set(gca,'FontSize',22)
 h_legend=legend('\alpha<1','\alpha=1')
 set(h_legend,'FontSize',22);  
 print('-depsc','graph/HUF.eps')
