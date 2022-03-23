% 
% Solve the Euler equation to determine the number of hires in a symmetric equilibrium
% Allows for possible layoffs
% np: past employment
% Rf: future recruiting costs
% a0: current technology
% hi: guess for number of hires
% 

function [h]=OBJEULER(np,a0,Rf,hi)

global delta eta c a alpha sigma_a omega rho_a s  B markup r varsigma  u_target th_target sigma w gamma
global q f u finv qinv uinv

u0=1-(1-s).*np;
w0=w.*a0^gamma;

RECR=@(hx)(hx>0).*a0.*c./q(finv(abs(hx)./u0)); % Recruiting costs

OBJF=@(hx)(RECR(hx)+w0-(1-s).*delta.*Rf-a0.*alpha.*(hx+(1-s).*np).^(alpha-1)).^2;
[res,val,exitflag]=fsolve(OBJF,hi,optimset('TolFun',10^(-13)));
h=res;
