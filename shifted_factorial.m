function result = shifted_factorial(c, j)

% Input:
%   c - the starting value of the shifted factorial
%   j - the number of terms in the product
%
% Output:
%   result - the value of (c)_j = c * (c+1) * ... * (c+j-1)
%

%   If j == 0, then (c)_j = 1 

    if j == 0
        result = 1; % Base case: (c)_0 = 1
    else
        result = 1; % Initialize the result
        for k = 0:j-1
            result = result * (c + k);
        end
    end
end