% 
% Shooting algorithm based on Fair and Taylor (1983)
% Generate solution to endogenous variable in current period and generate path for expectations of endogenous variables for period s on.
% Solve nonlinear RE system (based on certainty equivalence principle - see FT91)
% Y0: past variables
% s1: current state
% Gr: guess for expectations, includes a guess for current period (stochastic steady state)
% Y1: solution of the nl RE system -- endogenous variables in current period
% eII, eIII: final errors
% Er: final path of expectations from 2 to k+2h+1
% [Y1,er] is updated version of gr after 2 or 3 types of iteration (it may be longer if there are type III iterations)
% 

function [Y1,Er,kIII,eII,eIII,eII_record,eIII_record]=SHOOTING(Grt,Y0,s1,kIII)

global delta eta c a alpha sigma_a omega rho_a s  B r varsigma  u_target th_target sigma w gamma
global q f u finv qinv uinv
global ns ynum
global apos thpos npos mplpos hpos wpos Rpos  upos

%% Setup

tol2=5.*10^(-5); % Tolerance for convergence of algorithm (type II) - normal: 10^-6
tol3=10^(-3); % Tolerance for convergence of algorithm (type III) - normal: 10^-4
eIII_record=[];
eIII=1;
Thhat=50; % Must be big to force one iteration
Ehat=[Y0,Grt(:,1:kIII)]; % Current expectations - ynum*(kIII+1)
Yhat=[Y0,ones(ynum,kIII-1)]; % Current solutions - ynum*kIII

%% Iteration over length of path of expectations: type III iteration

while eIII>tol3
	eII_record=[];
	eII=1;

%% Iteration over path of expectations: type II iteration, forecast future endogenous variables

	while eII>tol2

%% Iteration over future time periods: type I iteration

		for i=2:kIII
			[Y]=SOLVESYS(Ehat(:,i-1),Ehat(:,i+1),Ehat(apos,i),Ehat(hpos,i));
			Yhat(:,i)=Y;
		end
		eII=max(max(abs(Yhat(:,1:end)-Ehat(:,1:end-1))));
		eII_record=[eII_record,eII];
		Ehat(:,1:end-1)=Yhat; % Update
	end

	figure(1)
	plot(Ehat(thpos,:))
	hold on
	eIII=abs(Yhat(thpos,2)-Thhat); % Compare estimates of theta
	Thhat=Yhat(thpos,2);
	eIII_record=[eIII_record,eIII];
	kIII=kIII+1; % Increase size of EP to check the result Y1 does not depend too much on final expectations
	Ehat=[Ehat,Grt(:,kIII)]; % Current expectations - ynum*(kIII+1) -- add one column
end

Y1=Yhat(:,2);
Er=Yhat(:,3:end); % Correct for additional entry
kIII=kIII-1;