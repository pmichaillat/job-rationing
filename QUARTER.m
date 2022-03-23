% % Transform the array from monthly frequency to quarterly frequency.
% Assume always starts in January
% Observation in columns
% 

function RES=QUARTER(DAT)

modulo=mod(size(DAT,1),3);
if modulo>0
 		  'number of quarter is not integer - suppress last months'
end
DAT=1/3*(DAT(1:3:end-2-modulo,:)+DAT(2:3:end-1-modulo,:)+DAT(3:3:end-modulo,:));
RES=DAT;
