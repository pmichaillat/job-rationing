%%==============================================================================
%% Create a time series for U froma  time series for TH
%%==============================================================================

function U_exact=CREATEU(THt,n0)

global s f
num=max(size(THt));
U_exact=[1-(1-s).*n0];

for i=1:num
	uold=U_exact(end);
	nold=(1-uold);
	nnew=nold+uold.*f(THt(i));
	unew=1-(1-s).*nnew;
 	U_exact=[U_exact,unew];
end
