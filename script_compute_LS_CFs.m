%% script_compute_LS_CFs  
%
% Description: 
% Script to compute LS-CFs, which are then investigated and compared
% in the other scripts 
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
dim = 1; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball, combi) 
weightFun = '1'; % weight function (1, C2k, sqrt)
F = 'trig'; % vecor space F (algebraic, trig, RBF) 
points = 'random'; % data points (equid, Halton, Sobol, Latin, random)  
d_max = 34; % max degree of exactness 

%% Start some loops 
for d = 0:d_max % start with d=1 for cubic PHS-RBF
            
    [d, d_max] % output to show the progress 
            
    %% Preliminary routines 
    omega = generate_weightFun( weightFun, dim); % set up weight function 
    [ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d ); % setup initial basis and moments 
    K = length(m); % number of basis functions 

    %% Compute a positive (nonnegative) LS-CF 
    [X, w, cond] = compute_LSCF( dim, domain, omega, basis, m, points ); % LS-CF 
    [X, w] = removeZeros( X, w ); % remove all zero weights and points

    %% Save the LS-CF 
    LS_CF = zeros(length(w),dim+2); 
    LS_CF(:,1:dim) = X; % data points 
    LS_CF(:,dim+1) = w; % cubature weights 
    LS_CF(1,dim+2) = K; % dimension of the function space 
    save( ['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat'], 'LS_CF' ); % safe matrix

end