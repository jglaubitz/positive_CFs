%% exp_highDim
% 
% Description: 
%  Function to compute the error of the CF applied to a high-dimensional test function on [-1,1]^dim
%
% INPUT: 
%  X :          	data points 
%  w :           	cubature weights 
%  CC :             number of experiments 
%
% OUTPUT:
%  error :  average errors for the Genz test functions 
%
% Author: Jan Glaubitz 
% Date: Jan 14, 2022 

function error = exp_highDim( X, w, CC )
    
	[N,dim] = size(X); % number of data points 
    error = 0; % error 
    
    for c=1:CC 

        %% set up Genz test functions 
        %a = ones(1,dim);
        a = rand(1,dim); b = rand(1,dim); % random parameters 
        fun = @(x) exp( -sum( (a.^2).*((x-b).^2) ) ); % Gaussian  
        I_ref = prod( sqrt(pi)./(2*a).*( erf( a.*b ) + erf( a.*(1-b) ) ) ); 
        fun_values = zeros(N,1); 
        for n=1:N 
            fun_values(n) = fun(X(n,:));
        end
        error = error + abs( I_ref - dot( w, fun_values ) ); % absolute error
             
    end
    
    % compute average errors 
    error = error/CC;
    
end