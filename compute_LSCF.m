%% compute_LSCFs  
%
% Description: 
%  Function to compute the LS-CF weights 
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 
% 
% INPUT: 
%  dim :        dimension 
%  domain :     (integration) domain 
%  omega :      weight function 
%  basis :      function space under consideration 
%  m :          vector of moments 
%  points :     type of data points 
%
% OUTPUT: 
%  X :          matrix that contains the data points 
%  w :          vector of cubature weights

function [ X, w] = compute_LSCF( dim, domain, omega, basis, m, points )

    K = length(m); % number of basis functions (dimension of F)

    %% routine to determine a nonnegative LS-CF 
    M = K; w_min = -1; 
    e = -1e-15; % tollerance to allow small rounding errors
    while w_min < e  
        
        %% data points, discrete weights, and coefficient matrix  
        Sample = generate_points( points, domain, dim, omega, M ); % data points and discrete weights
        R = sparse(diag(Sample.r)); % diagonal weight matrix 
        Phi = basis(Sample.coord); % initial Vandermonde matrix 
        
        %% Compute the LS weights 
        if rank(Phi) == K 
            A = Phi*sqrt(R); 
            v = lsqminnorm(A,m); % indirect computation using optimization tools 
            w = sqrt(R)*v;
            w_min = min(w); % their smallest value 
        end
        %M = 2*M; % increase the number of data points 
        M = M+K; % increase the number of data points
    end

    X = Sample.coord; % data points 
    
end