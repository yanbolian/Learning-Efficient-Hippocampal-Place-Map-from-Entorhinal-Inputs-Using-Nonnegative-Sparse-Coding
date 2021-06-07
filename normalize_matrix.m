function A_normalized = normalize_matrix(A, normalization_method, a_norm)
% This function normlizes the columns of A
% A: the matrix to be normalized
% normalization_method: 'L2 norm' or 'unit abs'; 
% a_norm: norm

% Set defaut value for 'normalizationMethod'
if ~exist('normalization_method', 'var')
    normalization_method = 'L2 norm';
end

% Set defaut value for 'aNorm'
if ~exist('a_norm', 'var')
    a_norm = 1;
end

% upper_bound = 1

% A small value to avoid zero division
epsilon = 1e-15;

% Normalize A
if isequal(normalization_method, 'L2 norm')
    % Normalize each column to l2-norm
    A_normalized = A * diag( a_norm ./ ( sqrt(sum(A.*A,1)) + epsilon ) );
elseif isequal(normalization_method,'unit abs')   
    % Normalize each column such that the sum of absolute values of column
    % elements equals to 'aNorm'
    A_normalized =  A * diag( a_norm ./ ( sum(abs(A),1)) + epsilon ); 
end

end