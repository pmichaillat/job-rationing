%
% Compute the ACF of all the time series in the first row of cell array DATA and place them in a cell array RES_ACF
% k_max is number of autocovariances to compute (so total number is k_max+1)
%

function RES_ACF=ACF(DATA,k_max)

numb=size(DATA,2);
RES_ACF=cell(5,numb);
for i=1:numb
	% Compute all the values of the ACF
	V=DATA{1,i};
	n=size(V,1);
	mu=mean(V,1);
	V2=V-mu;
	sigma=sum(V2.^2);
	n=min(k_max+1,n);
	Z=zeros(1,n);
	T=toeplitz(V2,Z);
	RES_ACF{1,i}=V2'*T./sigma;
	RES_ACF{2,i}=['ACF of ',DATA{2,i}];
	RES_ACF{3,i}='Lag k';
	RES_ACF{4,i}='ACF(k)';
	RES_ACF{5,i}=DATA{3,i}; % Identifiant	
end
end