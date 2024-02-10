%=======================================================================
%%Transform the array from weekly frequency to quarterly frequency.
%Assumes always starts in January
%=======================================================================

function RES=W2QUARTER(DAT)
%observation in rows

modulo=mod(size(DAT,1),12);
if modulo>0
 		  'number of quarter is not integer - please suppress last weeks'
end
DAT=1/3*(DAT(1:3:end-2-modulo,:)+DAT(2:3:end-1-modulo,:)+DAT(3:3:end-modulo,:));
RES=DAT;
