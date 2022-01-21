%% script_compute_Steinitz_CFs  
%
% Description: 
% Script to compute interpolatory CFs by Steinitz method
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
dim = 2; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball, combi) 
weightFun = '1'; % weight function (1, C2k, sqrt)
F = 'algebraic'; % vecor space F (algebraic, trig, cubic) 
points = 'Halton'; % data points (Halton, Sobol, Latin, random)  
d_max = 12; % maximum degree of exactness (usually 20 in 2d and 14 in 3d)

for d = 0:d_max  
            
    [d, d_max] % output to show the progress 
            
    %% Load the corresponding LS-CF 
    example = matfile(['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.LS_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights

    %% Basis, moments, and coefficient matrix 
    [ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d, points ); % setup basis and moments 
    Phi = basis(X); % coefficient matrix 
    K = length(m);

    %% Compute an interpolatory CF by Steinitz' method 
    tic; % start measuring time 
    X_Steinitz = X; w_Steinitz = w; % new variables 
    while K < length(w_Steinitz) 
        [dim, counter_points, d, K, length(w_Steinitz)]
        [X_Steinitz, w_Steinitz] = Steinitz( X_Steinitz, w_Steinitz, basis ); % apply Steinitz' Austauschsatz 
        [X_Steinitz, w_Steinitz] = removeZeros( X_Steinitz, w_Steinitz ); % remove all zero weights and points
    end 
    tEnd = toc; % stop measuring time 

    % Save the interpolatory CF by Steinitz' method 
    interpol_CF_Steinitz = zeros(length(w_Steinitz),dim+3); 
    interpol_CF_Steinitz(:,1:dim) = X_Steinitz; % data points 
    interpol_CF_Steinitz(:,dim+1) = w_Steinitz; % cubature weights 
    interpol_CF_Steinitz(:,dim+2) = N; % original number of data points 
    interpol_CF_Steinitz(1,dim+3) = tEnd; % time it took MATLAB to subsample 
    save( ['CFs/interpol_CF_Steinitz_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat'], 'interpol_CF_Steinitz' ); % safe matrix

end