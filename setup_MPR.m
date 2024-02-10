%%=================================================================
%% Setup the economic environment: real wage rigid, diminishing MPL, inelastic supply of labor; perfect competition; risk-neutral hh
%% Weekly frequency (1/12 of a quarter), 1/4 of a month
%%=================================================================

%clear all;close all;

global delta eta c a alpha sigma_a omega rho_a s z q f u B markup r varsigma finv qinv uinv u_target th_target w gamma ns

delta=(0.95).^(1./(12*4));      % discount factor
s=0.038./4;    % destruction rate
omega=0.933./4; 	%matching coeff
eta=0.5; 	%unemp-elast of matching function
c=.215;	 %recruiting costs -- .322 * wage
%steady-state
a=1; 		%LR technology

%------------------------------------------------------------------

%targets from JOLTS and CPS
u_target=0.058;% CPS -- 1964-2009

%epsi=9;% %product's elasticity
r=1-delta.*(1-s);
varsigma=s./(1-s);
q=@(x)omega.*x.^(-eta);
qinv=@(q)(q./omega).^(-1./eta);
finv=@(f)(f./omega).^(1./(1-eta));
f=@(x)omega.*x.^(1-eta);
uinv=@(ux)finv(varsigma*(1-ux)./ux);
u=@(th)varsigma./(varsigma+f(th));

%------------------------------------------------------------------------------
ft=omega.*(uinv(u_target)).^(1-eta);
alpha=1;
B=alpha.*(1-s).^(1-alpha);
gamma=.87;%real wage rigidity 
w=a./(1+r.*0.322./q(uinv(u_target)));
c=.322.*w;






