% 
% Solve the nonlinear system in period t
% Yp: past value of Y
% Y: solution of system ynum*1
% a0: current technology
% hi: initial guess for h
% Yf: expectation for Y in the next period
% 

function [Y]=SOLVESYS(Yp,Yf,a0,hi)

global delta eta c a alpha sigma_a omega rho_a s  B markup r varsigma  u_target th_target sigma w gamma
global q f u finv qinv uinv
global apos thpos npos mplpos hpos wpos Rpos ynum upos

Y=zeros(ynum,1);

np=Yp(npos);
u0=1-(1-s).*np;
w0=w.*a0^gamma;
Rf=Yf(Rpos);

%% Solve firm's Euler equation  

h0=OBJEULER(np,a0,Rf,hi); % h can be <0 if there are layoffs
	
%% Derive other solution values  

Y(apos)=a0;
Y(hpos)=h0;
Y(wpos)=w0;
Y(upos)=u0;
Y(thpos)=finv(max(h0,0)./u0);
Y(npos)=(1-s).*np+h0;
Y(mplpos)=a0.*alpha.*Y(npos).^(alpha-1);
Y(Rpos)=c.*a0./q(Y(thpos));