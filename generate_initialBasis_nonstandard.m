%% generate_initialBasis_nonstandard
%
% Description: 
%  Generates the initial basis and the corresponding moments for
%  polynomials on a nonstandard domain 
% 
% Author: Jan Glaubitz 
% Date: Jan 13, 2021 
%
% INPUT: 
%  d :          degree of exactness 
%
% OUTPUT: 
%  basis : vector-valued function with basis elements 
%  m :     corresponding moments 

function [ basis, m ] = generate_initialBasis_nonstandard( d )

    K = nchoosek(2 + d, 2); % binomial coefficient/ dimension
    alpha1 = zeros(K,1); % vector of exponents for x
    alpha2 = zeros(K,1); % vector of exponents for y
    m = zeros(K,1); % moments 
    
    %% exponents and basis 
    k = 1; 
   	for k1=0:1:d
   	for k2=0:1:d 
     	if k1+k2<=d 
         	alpha1(k) = k1; 
           	alpha2(k) = k2; 
           	k = k+1;
        end
    end
    end 
   	basis = @(x) x(:,1)'.^alpha1 .* x(:,2)'.^alpha2; % basis 
  
    
    %% moments of the monomials 
    fun_r = @(theta) 1 - (1/3)*sin(2*theta).^2;  
    k = 1; 
   	for k1=0:1:d
   	for k2=0:1:d 
       	if k1+k2<=d 
            % k1 and k2 need to be even
           	if mod(k1,2) == 0 && mod(k2,2) == 0 
                fun_moment = @(theta) (cos(theta).^k1).*(sin(theta).^k2).*(fun_r(theta).^(2+k1+k2)); 
                m(k) = integral(fun_moment,0,2*pi)/(2+k1+k2); 
            end
            k = k+1;
        end
    end
  	end
    
    
end