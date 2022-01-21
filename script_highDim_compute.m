%% script_highDim_compute
%
% Description: 
% Script to compute LS-CFs in high dimensions 
%
% Author: Jan Glaubitz 
% Date: Jan 13, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
dim = 10; % dimension (1,2,3,4) 
F = 'Gauss'; 
domain = 'cube'; % domain (cube, ball, combi) 
weightFun = '1'; % weight function (1, C2k, sqrt) 
points = 'random'; % data points (equid, Halton, Sobol, Latin, random)  
ep_ref = 0.01; 

I = logspace(0,3,10);
%% Start some loops 
for i = 1:10 % start with d=1 for cubic PHS-RBF
            
    [i, 10] % output to show the progress 
    d = floor(I(i));   
    
    %% Preliminary routines 
    ep = ep_ref*sqrt(d);
    omega = generate_weightFun( weightFun, dim); % set up weight function 
    [ basis, m ] = generate_initialBasis_Gauss_highDim( dim, domain, weightFun, d, ep, points  );
    K = length(m);
    
    %% Compute a positive LS RBF-CF 
    [X, w, cond] = compute_LSCF( dim, domain, omega, basis, m, points ); % LS-CF 
    N = length(w);
    % save 
    LS_CF = zeros(length(w),dim+2); 
    LS_CF(:,1:dim) = X; % data points 
    LS_CF(:,dim+1) = w; % cubature weights 
    LS_CF(1,dim+2) = K; % dimension of the function space 
    save( ['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat'], 'LS_CF' ); % safe matrix 
    
end