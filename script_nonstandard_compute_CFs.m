%% script_nonstandard_compute_CFs  
%
% Description: 
% Script to compute interpolatory and LS-CFs for a nonstandard domain  
%
% Author: Jan Glaubitz 
% Date: Jan 13, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
points = 'Halton'; % data points (equid, Halton, Sobol, Latin, random)  
d_max = 14; % max degree of exactness 

% Fixed parameters
dim = 2; % dimension 
F = 'algebraic'; % function space 
weightFun = '1'; % weight function 
domain = 'nonstandard'

%% Start some loops 
for d = 0:d_max % start with d=1 for cubic PHS-RBF
            
    [d, d_max] % output to show the progress 
            
    %% Preliminary routines 
    omega = generate_weightFun( weightFun, dim); % set up weight function 
    [ basis, m ] = generate_initialBasis_nonstandard( d );
    K = length(m); 
    
    %% Compute a positive LS-CF 
    [X, w, cond] = compute_LSCF( dim, domain, omega, basis, m, points ); % LS-CF 
    N = length(w);
    % save 
    LS_CF = zeros(length(w),dim+2); 
    LS_CF(:,1:dim) = X; % data points 
    LS_CF(:,dim+1) = w; % cubature weights 
    LS_CF(1,dim+2) = K; % dimension of the function space 
    save( ['CFs/LS_CF_',domain,'_d',num2str(d),'_',points,'.mat'], 'LS_CF' ); % safe matrix

    %% Compute and save an interpolatory CF by Steinitz' method 
    tic; % start measuring time 
    X_Steinitz = X; w_Steinitz = w; % new variables 
    while K < length(w_Steinitz) 
        [d, K, length(w_Steinitz)]
        [X_Steinitz, w_Steinitz] = Steinitz( X_Steinitz, w_Steinitz, basis ); % apply Steinitz' Austauschsatz 
        [X_Steinitz, w_Steinitz] = removeZeros( X_Steinitz, w_Steinitz ); % remove all zero weights and points
    end 
    tEnd = toc; % stop measuring time 
    % Save the interpolatory CF by Steinitz' method 
    interpol_CF = zeros(length(w_Steinitz),dim+3); 
    interpol_CF(:,1:dim) = X_Steinitz; % data points 
    interpol_CF(:,dim+1) = w_Steinitz; % cubature weights 
    interpol_CF(:,dim+2) = N; % original number of data points 
    interpol_CF(1,dim+3) = tEnd; % time it took MATLAB to subsample 
    save( ['CFs/interpol_CF_',domain,'_d',num2str(d),'_',points,'.mat'], 'interpol_CF' ); % safe matrix
    
end