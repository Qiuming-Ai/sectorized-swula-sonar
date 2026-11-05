function value = off_diagonal_max(matrix)
%OFF_DIAGONAL_MAX Maximum matrix entry excluding the main diagonal.
mask = ~eye(size(matrix,1),size(matrix,2));
values = matrix(mask);
if isempty(values), value = NaN; else, value = max(values); end
end
