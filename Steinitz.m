%% Steinitz 
%
% Description: 
%  Apply Steinitz' method to a positive CF to reduce the number of
%  data points 
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 
% 
% INPUT: 
%  X :      Matrix which contains the data points 
%  w :      Vector of cubature weights 
%  basis :  basis functions if the function space
%
% OUTPUT: 
%  X : Matrix which contains the data points 
%  v : Vector of new cubature weights (at least one is zero now)

function [ X, v ] = Steinitz( X, w, basis )

    %% Determine a suitable vector a 
    Phi = basis(X); % Vandermonde matrix 
   	Null = null(Phi); % Null space of Phi 
  	a = Null(:,1); % select the first basis element of the null space 
   	if max(a) <= 0 % no element is positive 
     	a = -a; % go over to -a 
    end
    
    %% Compute the new cubature weights 
    sigma = max( a./w ); % sigma 
    v = ( sigma*w - a )/sigma; % new cubature weights
    
end