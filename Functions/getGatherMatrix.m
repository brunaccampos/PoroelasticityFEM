function [L] = getGatherMatrix (rows,columns,e)
% 'Gather' matrix
% Input parameters: ndofs for rows and columns of gather matrix of element
% e

L = zeros(rows, columns);
for i = 1:rows
    for j = 1:columns
        if j == (rows-1)*(e-1)+i
            L(i,j) = 1;
        end
    end
end

end