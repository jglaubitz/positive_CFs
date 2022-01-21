%% script_compute_RBF_CFs  
%
% Description: 
% Script to compute interpolatory and LS RBF-CFs 
%
% Author: Jan Glaubitz 
% Date: Jan 12, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
dim = 2; % dimension (1,2,3) 
F = 'Gauss'; 
domain = 'cube'; % domain (cube, ball, combi) 
weightFun = '1'; % weight function (1, C2k, sqrt)
points = 'random'; % data points (equid, Halton, Sobol, Latin, random)  
K = 60; % max degree of exactness 
ep = 0.75; 

%% Start some loops 
for k = 1:K % start with d=1 for cubic PHS-RBF
            
    [k, K] % output to show the progress 
            
    %% Preliminary routines 
    omega = generate_weightFun( weightFun, dim); % set up weight function 
    %[ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, k, points ); % setup initial basis and moments 
    [ basis, m ] = generate_initialBasis_Gauss( dim, domain, weightFun, k, ep, points );
    
    %% Compute a positive LS RBF-CF 
    [X, w, cond] = compute_LSCF( dim, domain, omega, basis, m, points ); % LS-CF 
    N = length(w);
    % save 
    LS_CF = zeros(length(w),dim+2); 
    LS_CF(:,1:dim) = X; % data points 
    LS_CF(:,dim+1) = w; % cubature weights 
    LS_CF(1,dim+2) = K; % dimension of the function space 
    save( ['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(k),'_',points,'.mat'], 'LS_CF' ); % safe matrix

    %% Compute and save an interpolatory RBF-CF 
    % Compute 
    Sample = generate_points( points, domain, dim, omega, N ); % data points and discrete weights
    X = Sample.coord; 
    Phi = basis(X); 
    w = Phi\m; 
    % save 
    RBF_CF = zeros(length(w),dim+2); 
    RBF_CF(:,1:dim) = X; % data points 
    RBF_CF(:,dim+1) = w; % cubature weights 
    save( ['CFs/RBF_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(k),'_',points,'.mat'], 'RBF_CF' ); % safe matrix
    
end