%% generate_initialBasis_Gauss_highDim
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
%  ord:          number of centers 
%  ep :         shape parameter 
%
% OUTPUT: 
%  basis : vector-valued function with basis elements 
%  m :     corresponding moments 

function [ basis, m ] = generate_initialBasis_Gauss_highDim( dim, domain, weightFun, K, ep, points  )

    omega = generate_weightFun( weightFun, dim); % set up weight function
    Sample = generate_points( points, domain, dim, omega, K ); % data points and discrete weights
    X_RBF = Sample.coord; % data points
    m = zeros(K+1,1); % moments  
  
    %% 2D 
    if dim == 2
    
        % Compute the basis functions 
        rbf = @(r) exp(-(ep*r).^2); % thin plate spline (TPS)
        basis = @(x) [(x(:,1)').^0; rbf( sqrt( (x(:,1)'-X_RBF(:,1)).^2 + (x(:,2)'-X_RBF(:,2)).^2 ) )]; % RBF basis 
                
        %% Compute the moments 
        m(1) = 2^dim; 
        a = -1; b = 1; 
        for k=2:K 
            m1 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,1)) ) - erf( ep*(a-X_RBF(k,1)) ) ); % component in x direction
            m2 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,2)) ) - erf( ep*(a-X_RBF(k,2)) ) ); % component in x direction
            m(k) = m1*m2; % moments 
        end
             
    %% 3D
    elseif dim == 3
        
        % Compute the basis functions 
        rbf = @(r) exp(-(ep*r).^2); % thin plate spline (TPS)
        basis = @(x) [(x(:,1)').^0; rbf( sqrt( (x(:,1)'-X_RBF(:,1)).^2 + (x(:,2)'-X_RBF(:,2)).^2 + (x(:,3)'-X_RBF(:,3)).^2 ) )]; % RBF basis 
                
        %% Compute the moments 
        m(1) = 2^dim; 
        a = -1; b = 1; 
        for k=1:K 
            m1 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,1)) ) - erf( ep*(a-X_RBF(k,1)) ) ); % component in x direction
            m2 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,2)) ) - erf( ep*(a-X_RBF(k,2)) ) ); % component in x direction 
            m3 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,3)) ) - erf( ep*(a-X_RBF(k,3)) ) ); % component in x direction
            m(k+1) = m1*m2*m3; % moments 
        end
    
    %% 4D
    elseif dim == 4
        
        % Compute the basis functions 
        rbf = @(r) exp(-(ep*r).^2); % thin plate spline (TPS)
        basis = @(x) [(x(:,1)').^0; rbf( sqrt( (x(:,1)'-X_RBF(:,1)).^2 + (x(:,2)'-X_RBF(:,2)).^2 + (x(:,3)'-X_RBF(:,3)).^2  + (x(:,4)'-X_RBF(:,4)).^2 ) )]; % RBF basis ) )]; % RBF basis 
                
        %% Compute the moments 
        m(1) = 2^dim; 
        a = -1; b = 1; 
        for k=1:K 
            m1 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,1)) ) - erf( ep*(a-X_RBF(k,1)) ) ); % component in x direction
            m2 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,2)) ) - erf( ep*(a-X_RBF(k,2)) ) ); % component in x direction 
            m3 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,3)) ) - erf( ep*(a-X_RBF(k,3)) ) ); % component in x direction
            m4 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,4)) ) - erf( ep*(a-X_RBF(k,4)) ) ); % component in x direction
            m(k+1) = m1*m2*m3*m4; % moments 
        end
    
    %% 10D
    elseif dim == 10
        
        % Compute the basis functions 
        rbf = @(r) exp(-(ep*r).^2); % thin plate spline (TPS)
        basis = @(x) [(x(:,1)').^0; rbf( sqrt( (x(:,1)'-X_RBF(:,1)).^2 + (x(:,2)'-X_RBF(:,2)).^2 + ... 
            (x(:,3)'-X_RBF(:,3)).^2  + (x(:,4)'-X_RBF(:,4)).^2 + (x(:,5)'-X_RBF(:,5)).^2 + (x(:,6)'-X_RBF(:,6)).^2 + ... 
            (x(:,7)'-X_RBF(:,7)).^2  + (x(:,8)'-X_RBF(:,8)).^2 + (x(:,9)'-X_RBF(:,9)).^2 + (x(:,10)'-X_RBF(:,10)).^2 ... 
            ))]; % RBF basis ) )]; % RBF basis 
                
        %% Compute the moments 
        m(1) = 2^dim; 
        a = -1; b = 1; 
        for k=1:K 
            m1 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,1)) ) - erf( ep*(a-X_RBF(k,1)) ) ); % component in x1 direction
            m2 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,2)) ) - erf( ep*(a-X_RBF(k,2)) ) ); % component in x2 direction 
            m3 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,3)) ) - erf( ep*(a-X_RBF(k,3)) ) ); % component in x3 direction
            m4 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,4)) ) - erf( ep*(a-X_RBF(k,4)) ) ); % component in x4 direction
            m5 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,5)) ) - erf( ep*(a-X_RBF(k,5)) ) ); % component in x5 direction
            m6 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,6)) ) - erf( ep*(a-X_RBF(k,6)) ) ); % component in x6 direction 
            m7 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,7)) ) - erf( ep*(a-X_RBF(k,7)) ) ); % component in x7 direction
            m8 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,8)) ) - erf( ep*(a-X_RBF(k,8)) ) ); % component in x8 direction
            m9 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,9)) ) - erf( ep*(a-X_RBF(k,9)) ) ); % component in x9 direction
            m10 = (0.5/ep)*sqrt(pi)*( erf( ep*(b-X_RBF(k,10)) ) - erf( ep*(a-X_RBF(k,10)) ) ); % component in x10 direction 
            m(k+1) = m1*m2*m3*m4*m5*m6*m7*m8*m9*m10; % moments 
        end
        
    end
end