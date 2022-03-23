% 
% Construct Markov chain to replicate AR(1) process with parameters rho_m, conditional variance sigma_e, and mean 0
% ns states
% Based on the algorithm described in Adda and Cooper (MIT Press, 2002)
% 

function [Z,PI,EPS]=MAKEMC(ns,rho_m,sigma_e)

prec=1000;
sigma_m=(sigma_e.^2./(1-rho_m^2))^(0.5);
EPS=[-25,sigma_m.*norminv(([2:ns]-1)./ns,0,1),25];
Z=-ns.*sigma_m.*normpdf(EPS(2)./sigma_m,0,1);
Z=[Z,ns.*sigma_m.*(normpdf(EPS(2:ns-1)./sigma_m,0,1)-normpdf(EPS(3:ns)./sigma_m,0,1))];
Z=[Z,ns.*sigma_m.*normpdf(EPS(ns)./sigma_m,0,1)];

for i=1:ns
   if (i==1) | (i==ns)
      pas=(EPS(i+1)-EPS(i))./(200*prec);
    else
      pas=(EPS(i+1)-EPS(i))./prec;
    end
  xs=[EPS(i):pas:EPS(i+1)];
  for j=1:ns
    ys=exp(-xs.^2./(2.*sigma_m.^2)).*(normcdf((EPS(j+1)-rho_m.*xs)./sigma_e,0,1)-normcdf((EPS(j)-rho_m.*xs)./sigma_e,0,1));
    PI(i,j)=ns./(2.*pi.*sigma_m.^2).^(1/2).*trapz(xs,ys);
  end
end
