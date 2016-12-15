function [delZ] = abs2inc(Z)
% changes data matrix Z from absolute value to incremental
% for incremental model identification

[n_sam n_var] = size(Z);	% obtain dimensions of original matrix

delZ = zeros(n_sam, n_var);	% define output matrix size (same as input)

delZ(2:end,:) = Z(2:end,:) - Z(1:end-1,:);

end