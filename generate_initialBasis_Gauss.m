%% generate_initialBasis_Gauss
%
% Description: 
%  Function that generates the basis and the corresponding moments when F is the
%  space of quintic-RBFs corresponding to a set of centers X_RBF
%
% Author: Jan Glaubitz 
% Date: Jan 12, 2021 
%
% INPUT: 
%  dim :        dimension 
%  domain :     domain 
%  weightFun :  weight function 
%  K :          number of centers 
%  points :     type of data points 
%
% OUTPUT: 
%  basis : vector-valued function with basis elements 
%  m :     corresponding moments 

function [ basis, m ] = generate_initialBasis_Gauss( dim, domain, weightFun, K, ep, points  )

    omega = generate_weightFun( weightFun, dim); % set up weight function
    Sample = generate_points( points, domain, dim, omega, K ); % data points and discrete weights
    X_RBF = Sample.coord; % data points
    m = zeros(K+1,1); % moments 
    
    %% exponents and basis 
    if dim == 2 
        if strcmp( domain, 'cube') 
            if strcmp( weightFun, '1') 
                
                %% Compute the basis functions 
                rbf = @(r) exp(-(ep*r).^2); % thin plate spline (TPS)
                %basis = @(x) rbf( sqrt( (x(:,1)'-X_RBF(:,1)).^2 + (x(:,2)'-X_RBF(:,2)).^2 ) ); % RBF basis 
                basis = @(x) [(x(:,1)').^0; rbf( sqrt( (x(:,1)'-X_RBF(:,1)).^2 + (x(:,2)'-X_RBF(:,2)).^2 ) )]; % RBF basis 
                
                %% Compute the moments 
                m(1) = 4; 
                a = -1; b = 1; 
                for k=1:K 
                    mx = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,1)) ) - erf( ep*(a-X_RBF(k,1)) ) ); % component in x direction
                    my = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,2)) ) - erf( ep*(a-X_RBF(k,2)) ) ); % component in x direction
                    %m(k) = mx*my; % moments 
                    m(k+1) = mx*my; % moments 
                end
            else 
                error('Desired weight function not yet implemented!')
            end
        else 
           error('Desired domain not yet implemented!') 
        end
    else
        error('Desired dimension not yet implemented!') 
    end
        
end