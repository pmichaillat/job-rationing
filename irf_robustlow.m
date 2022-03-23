% 
% Compute IRFs of log-linearized DSGE model with job rationing and matching frictions
% Assume low recruiting cost
% 

clear all;close all;
setupsimul_robustlow;
T=400; % Length of IRF

%% Compute IRFs

shock=zeros(xeq,1);
shock(a_pos)=-sigma_a;
imptech=BMAT*shock;
for i=2:T
    imptech=[imptech AMAT*imptech(:,i-1)]; % Use reduced form solution: x_t = amat * x_{t-1} + b * shock
end

% Include frictional and rationing unemployment
StS          = SS(1:xvari);
ucpos        = xvari+1;ufpos=xvari+2;
StS(ucpos,1) = UC_bar;
StS(ufpos,1) = UF_bar;

IRF = imptech(1:xvari,:); % Store impulse responses
IRF = [IRF;-IRF(apos,:).*(1-gamma)./(1-alpha).*(1-UC_bar)./UC_bar];
IRF = [IRF;U_bar./UF_bar.*IRF(upos,:)-UC_bar./UF_bar.*IRF(end,:)];

%% Plot IRFs

top=250;
rang=[1:top];
gra=[apos,ypos,wpos,hpos,thpos,upos,9,10]; 
nam={'Technology','Output','Wage','Number of hires','Vacancy-unemployment ratio','Unemployment','Rationing unemployment','Frictional unemployment'};
ba=[1,Y_bar,W_bar,H_bar,TH_bar,U_bar,UC_bar,UF_bar];
la=[-0.0035,-0.0035,-0.0035,-0.025,-0.06,0,0,-0.15];
ma=[0,0,0,0.001,0,0.025,2.5,0];

figure(1)
clf
for j=1:8
  subplot(4,2,j)
  hold on
  plot(rang,IRF(gra(j),rang),'-b','LineWidth',3)
  set(gca,'YGrid','on','XTickLabel','','FontSize',14)
  xlim([0 top]) ;
  ylim([la(j) ma(j)]);
  title(nam{j},'FontSize',14)
end
print('-depsc','IRF_robustlow.eps')