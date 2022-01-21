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
% Date: Jan 21, 2021 
% 

function [ basis, m ] = generate_initialBasis_trigPoly( dim, domain, weightFun, d )
    
    % vectors of exponents
    alpha = (1:d)'; % vector of exponents for x 
    K = 1 + 2*d; 
    
    %% basis 
    if dim == 1 
        basis = @(x) [(x').^0; sin( alpha*pi*x' ); cos( alpha*pi*x' ) ];
    else 
    	error('Desired dimension not yet implemented!') 
    end
    
    %% moments of the monomials 
    m = zeros(K,1); % moments
    % cube 
    if strcmp( domain, 'cube') & strcmp( weightFun, '1') 
        m(1) = 2; 
            
    % else 
    else
        error('Desired domain or weight function not yet implemented!')
    end
    
end