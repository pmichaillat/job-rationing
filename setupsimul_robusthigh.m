% 
% Set up for simulation of log-linear model with job rationing for IRF and simulated moments
% Assume high recruiting cost
% 

global  W_bar C_bar TH_bar  N_bar H_bar U_bar UC_bar UF_bar Y_bar 
global rho_a sigma_a nsample
global wpos cpos thpos npos hpos upos ypos apos  a_pos

%% Calibration

setup_robusthigh; % Model parameters
[W_bar,C_bar,TH_bar,N_bar,H_bar,U_bar,UC_bar,UF_bar,Y_bar]=STEADYLL(w,gamma); % Steady state
nsample=182;
[rho_a,sigma_a,ax]=TECHNO(alpha,nsample); % Technology process

% Positions of variables
wpos	=2;					% wage
cpos   =  3;               % consumption
thpos   =  1;               % market tightness
npos   =  4;               % employment
hpos   =  5;               % hiring rate
upos   =  6;                % unemployment rate
ypos  =  7;               % output
apos   =  8;               % technology 
a_pos   = 9;               % technology shock

%% Solve log-linear model

[AMAT,BMAT,SS,xeq,xvari]=LIN_DSGE(w,gamma,alpha); % To compute and plot impulse response use reduced form solution: x_t = amat * x_{t-1} + b * shock