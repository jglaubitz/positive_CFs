%% generate_initialBasis
%
% Description: 
%  Generates the initial basis and the corresponding moments 
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 
% 
% INPUT: 
%  dim :        dimension 
%  domain :     domain 
%  weightFun :  weight function 
%  F :          vector space F 
%  d :          degree of exactness 
%  points :     type of data points
%
% OUTPUT: 
%  basis : vector-valued function with basis elements 
%  m :     corresponding moments 

function [ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d, points )

    K = nchoosek(dim + d, dim); % dimension of space of polynomials of total degree up to d

    if strcmp( F, 'algebraic')  
        [ basis, m ] = generate_initialBasis_algPoly( dim, domain, weightFun, d ); 
    elseif strcmp( F, 'trig')
        [ basis, m ] = generate_initialBasis_trigPoly( dim, domain, weightFun, d );
    elseif strcmp( F, 'cubic') 
        [ basis, m ] = generate_initialBasis_cubic( dim, domain, weightFun, K, points );
    else
        error('Desired vector space F not yet implemented!')
    end
    
end