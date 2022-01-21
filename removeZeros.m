%% removeZeros 
%
% Description: 
%  Function to remove zero weights and the corresponding data points from a
%  CF 
% 
% INPUT: 
%  X : Matrix which contains the data points 
%  w : Vector of cubature weights
%
% OUTPUT: 
%  Y : Matrix which contains the data points 
%  v : Vector of cubature weights (no zeros anymore)
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 

function [ Y, v] = removeZeros( X, w )

    [N,dim] = size(X); 
    n = 1; 
    Y = X; v = w; % start from the original CF
    
    e = 1e-15; % tollerance to allow small rounding errors
    while n <= N 
        if abs( v(n) ) <= e % zero weight 
            Y(n,:) = []; % remove the data points
            v(n) = []; % remove the weight 
            N = N-1; 
        else 
            n = n+1; 
        end
    end 
    
end