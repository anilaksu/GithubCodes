function [lu, pvt] = LUfactor ( A )

%LUFACTOR     compute the LU decomposition for the matrix A
%
%     calling sequence:
%             [lu, pvt] = LUfactor ( A )
%
%     input:
%             A       coefficient matrix for linear system
%                     (matrix must be square)
%
%     outputs:
%             lu      matrix containing LU decomposition of input matrix
%                     A - the L matrix in the decomposition consists of
%                     1's along the main diagonal together with the 
%                     strictly lower triangular portion of the matrix lu;
%                     the U matrix in the decomposition is the upper
%                     triangular portion of the matrix lu
%             pvt     vector which indicates the permutation of the rows
%                     performed during the decomposition process
%
%     NOTE:
%             this routine performs partial pivoting during the 
%             decomposition process
%
%             the system Ax = b can be solved by first applying LUfactor
%             to the coefficient matrix A and then applying the companion
%             routine, LUsolve, for each right-hand side vector b
%

[nrow ncol] = size ( A );
if ( nrow ~= ncol )
   disp ( 'LUfactor error: Square coefficient matrix required' );
   return;
end;

%
%   initialize row pointers
%

for i=1:nrow
    pvt(i) = i;
end;

for i = 1 : nrow - 1

%
%   partial pivoting
%

% This part was removed by Pete Diamessis

%
%   terminate if matrix is singular
%

	if ( A(pvt(i),i) == 0 ) 
	   disp ( 'LUfactor error: coefficient matrix is singular' );
	   lu = A;
	   return
	end;
	
%
%   elimination steps
%

    for j = i+1 : nrow
	    m = -A(pvt(j),i) / A(pvt(i),i);
		A(pvt(j),i) = -m;
		A(pvt(j), i+1:nrow) = A(pvt(j), i+1:nrow) + m * A(pvt(i), i+1:nrow);
    end;
end;

lu = A;

