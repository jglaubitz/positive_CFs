%% generate_initialBasis_trigPoly
%
% Description: 
%  Function to generate initial basis and the corresponding moments when F is the
% space of trigonometric polynomials of degree at most d
% 
% INPUT: 
%  dim :        dimension 
%  domain :     domain 
%  weightFun :  weight function 
%  d :          degree of exactness 
%
% OUTPUT: 
%  basis : vector-valued function with basis elements 
%  m :     corresponding moments 
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 
% 

function [ basis, m ] = generate_initialBasis_trigPoly( dim, domain, weightFun, d )

    % dimension of the vector space 
    if d == 0 
        K = 1; 
    elseif d == 1 
        K = 5;
    else
        K = 0; 
        for j=0:dim 
            K = K + (2^j)*nchoosek(dim,j)*nchoosek(d,j);
        end 
    end
    
    % vectors of exponents
    alpha = zeros(K,1); % vector of exponents for x
    beta = zeros(K,1); % vector of exponents for y
    gamma = zeros(K,1); % vector of exponents for z  
    
    %% exponents 
    if dim == 1 
        alpha = (-d:d)'; 
    elseif dim == 2 
        k = 1; 
        for k1=-d:d
        for k2=-d:d 
            if abs(k1)+abs(k2)<=d 
                alpha(k) = k1; 
                beta(k) = k2; 
                k = k+1;
            end
        end
        end 
    elseif dim == 3 
    	k = 1; 
        for k1=-d:d
        for k2=-d:d 
        for k3=-d:d
        	if abs(k1)+abs(k2)+abs(k3)<=d 
            	alpha(k) = k1; 
                beta(k) = k2;
                gamma(k) = k3; 
                k = k+1;
            end
        end
        end
        end 
    else 
    	error('Desired dimension not yet implemented!') 
    end   
    
    %% delete negative counterparts 
    k = 1; 
    while k <= K 
        l = k+1; 
        while l <= K 
            if alpha(k) == -alpha(l) && beta(k) == -beta(l) && gamma(k) == -gamma(l)
                alpha(l) = []; 
                beta(l) = []; 
                gamma(l) = []; 
                K = K-1; 
                l = K+1;
            end 
            l = l+1;
        end
        k = k+1;
    end
    
    %% basis 
    if dim == 1  
        basis = @(x) exp( 2*pi*1i*alpha*(x') ) + exp( -2*pi*1i*alpha*(x') ); % basis 
    elseif dim == 2 
        basis = @(x) exp( 2*pi*1i*alpha*(x(:,1)') ).*exp( 2*pi*1i*beta*(x(:,2)') )... 
                + exp( -2*pi*1i*alpha*(x(:,1)') ).*exp( -2*pi*1i*beta*(x(:,2)') ); % basis
    elseif dim == 3 
        basis = @(x) exp( 2*pi*1i*alpha*(x(:,1)') )... 
                    .*exp( 2*pi*1i*beta*(x(:,2)') )... 
                    .*exp( 2*pi*1i*gamma*(x(:,3)') ) + ... 
                     exp( -2*pi*1i*alpha*(x(:,1)') )... 
                    .*exp( -2*pi*1i*beta*(x(:,2)') )... 
                    .*exp( -2*pi*1i*gamma*(x(:,3)') ); % basis
    else 
    	error('Desired dimension not yet implemented!') 
    end
    
    %% moments of the monomials 
    m = zeros(K,1); % moments
    % cube 
    if strcmp( domain, 'cube') 
        % omega = 1
        if strcmp( weightFun, '1') 
            if dim <= 3  
                for k=1:K 
                    if alpha(k) == 0 && beta(k) == 0 && gamma(k) == 0
                        m(k) = 2*(2^dim); % moment 
                    end
                end
            else 
                error('Desired dimension not yet implemented!') 
            end
        else 
            error('Desired weight function not yet implemented!')
        end
            
    % else 
    else
        error('Desired domain not yet implemented!')
    end
    
end