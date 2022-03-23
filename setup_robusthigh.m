% 
% Setup the economic environment: real wage rigid, diminishing marginal returns to labor, inelastic supply of labor; perfect competition; risk-neutral households
% Weekly frequency
% Assume high recruiting cost
% 

global delta eta c a alpha omega s  B r varsigma  u_target th_target w gamma ar
global q f u finv qinv uinv ur
global ns

%% Estimated parameters

delta=(0.95).^(1./(12*4)); % Discount factor
s=0.038./4; % Job-destruction rate
omega=0.933./4; % Matching efficacy
eta=0.5; % Unemployment elasticity of matching function
a=1; % Long-run technology
u_target=0.058; % CPS, 1964–2009
r=1-delta.*(1-s);
varsigma=s./(1-s);
gamma=.7; % Real-wage rigidity

%% Functions

% Use Cobb-Douglas matching function
q=@(x)omega.*x.^(-eta);
qinv=@(q)(q./omega).^(-1./eta);
finv=@(f)(f./omega).^(1./(1-eta));
f=@(x)omega.*x.^(1-eta);
uinv=@(ux)finv(varsigma*(1-ux)./ux);
u=@(th)varsigma./(varsigma+f(th));
th_target=uinv(u_target);
n_target=(1-u_target)./(1-s);

%% Calibration 

ls=0.66;
ccoeff=0.54;
alpha=ls.*(1+r.*ccoeff./q(th_target));
B=alpha.*(1-s).^(1-alpha); % Useful parameter
w=ls./n_target^(1-alpha);
c=ccoeff.*w;
ar=(w./alpha).^(1./(1-gamma));
ur=@(ax) max(1-((alpha./w).^(1./(1-alpha)).*ax.^((1-gamma)./(1-alpha))),0);