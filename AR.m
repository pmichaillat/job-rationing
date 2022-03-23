%	
% TS is a time series (column)
% Conduct AR(p)
%

function [Rho,std]=AR(TS,p)

T=max(size(TS)); % Use all the time series available. Condition on first p observations.
M=toeplitz(TS);
s2=0;
E_hat=[];

% Compute linear regressions on j lags for all the lags
% Output an (unbiased) estimate of the autoregression on p lags

X=M(1:p,p:end-1)';XX_inv=(X'*X)^(-1);
E_hat=(eye(T-p)-X*XX_inv*X')*TS(p+1:end);	
std=sqrt(E_hat'*E_hat./(T-p));
Rho=XX_inv*X'*TS(p+1:end);

