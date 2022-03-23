% 
% Setup the economic environment: real wage rigid, diminishing marginal returns to labor, inelastic supply of labor; perfect competition; risk-neutral households
% Weekly frequency
% 

clear all;close all;

global delta eta c a alpha sigma_a omega rho_a s z q f u B markup r varsigma finv qinv uinv u_target th_target w gamma ns

delta=(0.95).^(1./(12*4)); % Discount factor
s=0.038./4; % Job-destruction rate
omega=0.933./4; % Matching efficacy
eta=0.5; % Unemployment elasticity of matching function
c=.215; % Recruiting cost: 0.322 * wage
a=1; % Long-run technology
u_target=0.058; % CPS, 1964–2009
r=1-delta.*(1-s);
varsigma=s./(1-s);
q=@(x)omega.*x.^(-eta);
qinv=@(q)(q./omega).^(-1./eta);
finv=@(f)(f./omega).^(1./(1-eta));
f=@(x)omega.*x.^(1-eta);
uinv=@(ux)finv(varsigma*(1-ux)./ux);
u=@(th)varsigma./(varsigma+f(th));
ft=omega.*(uinv(u_target)).^(1-eta);
alpha=0.666;
B=alpha.*(1-s).^(1-alpha);
gamma=.7; % Real-wage rigidity 
w=0.671;	

% Parameters for model with rigid wage
wMPS=a./(1+r.*0.322./q(uinv(u_target)));
cMPS=.322.*wMPS;

% Parameters for canonical model
th_target=uinv(u_target);
OBJ=@(betax)(cMP/q(th_target).*r+cMP.*betax.*(1-s).*delta.*th_target-(1-betax).*a).^2;
[res,val,exitflag]=fsolve(OBJ,.5,optimset('TolFun',10^(-16)));
betaMP=res;
cMP=cMPS;

% Parameters for model with diminishing returns
ls=0.66; % Standard labor share
nbar=(1-u_target)./(1-s);
kappa=(r.*0.322./q(th_target)+1).*ls;
betaSZ=ls./(kappa+delta.*(1-s).*th_target.*0.322.*ls);
alphaSZ=(kappa-betaSZ.*kappa)./(1-betaSZ.*kappa);
wbar=ls.*nbar.^(alpha-1);
cSZ=0.322.*wbar;