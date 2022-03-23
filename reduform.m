%
% [nonsing,b] = reduform(q,qrows,qcols,bcols,neq,b,condn);
%
% Compute reduced-form coefficient matrix, b.
%

function [nonsing,b] = reduform(q,qrows,qcols,bcols,neq,condn);

left = 1:qcols-qrows;
right = qcols-qrows+1:qcols;

nonsing = rcond(q(:,right)) > condn;

if(nonsing)
   q(:,left) = -q(:,right)\q(:,left);
end

b = q(1:neq,1:bcols);

return
