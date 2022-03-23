% 
% Compute steady state of decentralized economy
% solve for steady-state of  given w, gamma
% w: rigid, steady-state wage level
% gamma: real wage rigidy
% 

function [Y]=STEADYGE(w,gamma)

global delta eta c a alpha beta sigma_a omega s rho_a z q f u markup sigma
global apos thpos npos mplpos hpos wpos Rpos ynum upos

Y=zeros(ynum,1);

% Theta in steady-state in both cases 

TH=FINDTH(w,gamma);

% Derive other steady-state values 

Y(apos)=a;
Y(thpos)=TH;
Y(upos)=u(TH);
Y(npos)=(1-u(TH))./(1-s);
Y(hpos)=s*Y(npos);
Y(Rpos)=c*a./q(TH);
Y(wpos)=w.*a^gamma;
Y(mplpos)=alpha.*a.*Y(npos).^(alpha-1);