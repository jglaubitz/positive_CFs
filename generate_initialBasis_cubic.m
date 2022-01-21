%% generate_initialBasis_cubic
%
% Description: 
%  Function that generates the basis and the corresponding moments when F is the
%  space of cubic-RBFs corresponding to a set of centers X_RBF
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 
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

function [ basis, m ] = generate_initialBasis_cubic( dim, domain, weightFun, K, points  )

    omega = generate_weightFun( weightFun, dim); % set up weight function
    Sample = generate_points( points, domain, dim, omega, K ); % data points and discrete weights
    X_RBF = Sample.coord; % data points
    m = zeros(K,1); % moments 
    
    %% exponents and basis 
    if dim == 2 
        if strcmp( domain, 'cube') 
            if strcmp( weightFun, '1') 
                
                %% Compute the basis functions 
                rbf_cubic = @(r) r.^3; % thin plate spline (TPS)
                basis = @(x) rbf_cubic( sqrt( (x(:,1)'-X_RBF(:,1)).^2 + (x(:,2)'-X_RBF(:,2)).^2 ) ); % RBF basis 
                
                %% Compute the moments 
                I_tr = @(u,v) (u/40)*( ... 
                    3*u^4*asinh(v/u) + ... 
                    v*( 5*u^2 + 2*v^2 )*sqrt( u^2 + v^2 ) ...
                ); % reference integral
                % compute moments 
                for k=1:K 
                    % shifted edges of the rectangle 
                    a = -1; b = 1; c = a; d = b; % we assume the domain [-1,1]^2
                    a_tilde = abs( a - X_RBF(k,1) ); 
                    b_tilde = abs( b - X_RBF(k,1) ); 
                    c_tilde = abs( c - X_RBF(k,2) ); 
                    d_tilde = abs( d - X_RBF(k,2) ); 
                
                    % partition rectangle in 8 right triangles and compute the
                    % corresponding integrals 
                    I(1) = I_tr(b_tilde,d_tilde); 
                    I(2) = I_tr(d_tilde,b_tilde); 
                    I(3) = I_tr(d_tilde,a_tilde); 
                    I(4) = I_tr(a_tilde,d_tilde); 
                    I(5) = I_tr(a_tilde,c_tilde); 
                    I(6) = I_tr(c_tilde,a_tilde); 
                    I(7) = I_tr(c_tilde,b_tilde); 
                    I(8) = I_tr(b_tilde,c_tilde); 
                    I(isnan(I))=0; % set all NaN values to zero;
                    % sum these up to get the moment 
                    m(k) = ( b_tilde*d_tilde ~= 0 )*(I(1)+I(2)) + ... 
                           ( a_tilde*d_tilde ~= 0 )*(I(3)+I(4)) + ...
                           ( a_tilde*c_tilde ~= 0 )*(I(5)+I(6)) + ...
                           ( b_tilde*c_tilde ~= 0 )*(I(7)+I(8));
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