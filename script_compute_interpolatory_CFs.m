%% script_compute_interpolatory_CFs  
%
% Description: 
% Script to compute interpolatory CFs from a given LS-CF by different
% subsampling methods 
%
% Note that there is a separate script for Steinitz, which is used for most
% of the tests, while basis pursuit and NNLS do not always work. 
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
d = 0; % degree of exactness (0,1,2,...)
points = 'Halton'; % data points (Halton, Sobol, random)  

%% Start some loops 
for d = 0:14  
            
	%[d, 14] % output to show the progress 
            
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
        [d, K, length(w_Steinitz)]
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

    %% Compute an interpolatory CF by basis pursuit (linear programming) 
    tic; % start measuring time 
    X_BP = X; w_BP = w; % new variables 
    options = optimoptions('linprog','Display','none');
    w_BP = linprog( ones(length(w),1), [], [], Phi, m, zeros(length(w),1), [], options ); % BP weights 
    [X_BP, w_BP] = removeZeros( X_BP, w_BP ); % remove all zero weights and points
    tEnd = toc; % stop measuring time 

    % Save the interpolatory CF by basis pursuit (linear programming) 
    interpol_CF_BP = zeros(length(w_BP),dim+3); 
    interpol_CF_BP(:,1:dim) = X_BP; % data points 
    interpol_CF_BP(:,dim+1) = w_BP; % cubature weights 
    interpol_CF_BP(:,dim+2) = N; % original number of data points 
    interpol_CF_BP(1,dim+3) = tEnd; % time it took MATLAB to subsample 
    save( ['CFs/interpol_CF_BP_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat'], 'interpol_CF_BP' ); % safe matrix

    %% Compute an interpolatory CF by nonnegative least-squares (NNLS)
    tic; % start measuring time 
    X_NNLS = X; w_NNLS = w; % new variables 
    options = optimset('Display', 'none');
    w_NNLS = lsqnonneg( Phi, m, options ); % NNLS weights 
    [X_NNLS, w_NNLS] = removeZeros( X_NNLS, w_NNLS ); % remove all zero weights and points
    tEnd = toc; % stop measuring time 

    % Save the interpolatory CF by nonnegative least-squares (NNLS)
    interpol_CF_NNLS = zeros(length(w_NNLS),dim+3); 
    interpol_CF_NNLS(:,1:dim) = X_NNLS; % data points 
    interpol_CF_NNLS(:,dim+1) = w_NNLS; % cubature weights 
    interpol_CF_NNLS(:,dim+2) = N; % original number of data points 
    interpol_CF_NNLS(1,dim+3) = tEnd; % time it took MATLAB to subsample 
    save( ['CFs/interpol_CF_NNLS_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat'], 'interpol_CF_NNLS' ); % safe matrix

end