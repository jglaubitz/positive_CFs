%% script_LS_MC_
%
% Description: 
% Script for LS-MC formulas in high dimensions 
%
% Author: Jan Glaubitz 
% Date: Jan 17, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
dim = 10; % dimension (1,2,3,4) 
F = 'algebraic'; % function space 
domain = 'cube'; % domain (cube, ball, combi) 
weightFun = '1'; % weight function (1, C2k, sqrt) 
points = 'Halton'; % data points (equid, Halton, Sobol, Latin, random)  

%% Prepare script 

% Test function and its integral 
CC = 40; % number of tests for Genz 
noise_level = 0; % noise added to the function values 

% initiate vectors 
error_MC = []; error_LS1 = []; error_LS2 = []; % errors

NN = logspace(2,5,10)';
% Start some loops 
for n = 1:length(NN) % start with d=1 for cubic PHS-RBF
            
    [n, length(NN)] % output to show the progress 
    N = floor(NN(n));   
    
    %% QMC without correction 
    omega = generate_weightFun( weightFun, dim); % set up weight function 
    Sample = generate_points( points, domain, dim, omega, N );
    X = (Sample.coord+1)/2;
    w = ones(N,1)/N;
    error = exp_highDim( X, w, CC ); % average errors   
    error_MC = [error_MC; error]; % errors corresponding to Genz' test functions 
    
    %% QMC with 1st-order correction (m=1) 
    d=1; % inlcude linear polynomials 
    [ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d, points ); % compute basis 
    Phi = basis(Sample.coord); % Vandermonde matrix 
    w = lsqminnorm(Phi,m); % LS weights on the cube [-1,1]^dim
    w = w/(2^dim); % LS weights on the cube [0,1]^dim 
    error = exp_highDim( X, w, CC ); % average errors   
    error_LS1 = [error_LS1; error]; % errors corresponding to Genz' test functions 
    
    %% QMC with 1st-order correction (m=2) 
    d=2; % inlcude linear polynomials 
    [ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d, points ); % compute basis 
    Phi = basis(Sample.coord); % Vandermonde matrix
    w = lsqminnorm(Phi,m); % LS weights on the cube [-1,1]^dim
    w = w/(2^dim); % LS weights on the cube [0,1]^dim 
    error = exp_highDim( X, w, CC ); % average errors   
    error_LS2 = [error_LS2; error]; % errors corresponding to Genz' test functions 
    
end

%% Plot figures for compparison - N vs error
figure(1) 
p = plot( NN,error_MC,'k^', NN,error_LS1,'rs', NN,error_LS2,'bo' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 3, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[f] - C_N[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('QMC','LS-QMC, $m=1$','LS-QMC, $m=2$','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
