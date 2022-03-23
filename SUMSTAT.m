% 
% Compute summary stat of matrix X (each row is an observation and each column is a variable.)
% 


function [moy,dev,autoc,M]=SUMSTAT(X)
	
moy=mean(X,1);
dev=std(X,1);
autoc=[];
for j=1:size(X,2)
autoc=[autoc,AUTOCORREL(X(:,j),2)];
end
autoc=autoc(2,:); % Keep serial correlation
M=corrcoef(X);
