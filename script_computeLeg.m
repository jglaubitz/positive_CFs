%% script_compute_Steinitz_CFs  
%
% Description: 
% Script to compute (transformed) product Legendre rules
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
dim = 2; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball) 
d_max = 20; % max total degree

for d=0:d_max
    
    % load NNI-CF to use approximately the same number of data points 
    example = matfile(['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_Halton.mat']);
    C = example.LS_CF; 
    [ M, aux] = size(C);
    n = ceil(M^(1/dim));
    
    % compute the Legendre rule 
    [X, w_Leg, d_Leg, K_Leg ] = compute_LegendreRule( dim, domain, n );
    
    % save points, weights, and d and K in a matrix
    Leg_CF = zeros(n^dim,dim+3); % initiate 
    Leg_CF(:,1:dim) = X; % store data points
    Leg_CF(:,dim+1) = w_Leg; % store cubature weights 
    Leg_CF(1,dim+2) = d_Leg; % store d
    Leg_CF(1,dim+3) = K_Leg; % store K      
    save( ['CFs/Leg_CF_dim',num2str(dim),'_',domain,'_d',num2str(d),'.mat'], 'Leg_CF' ); % safe matrix   
    
end