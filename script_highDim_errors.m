%% script_highDim_errors
%
% Description: 
%  Script to compare the accuracy and stbility of interpolatory RBF-CFs and
%  LS RBF-CFs 
%  
% Author: Jan Glaubitz 
% Date: Jan 13, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters 
dim = 10; % dimension 
domain = 'cube'; % domain (cube, ball, combi) 
weightFun = '1'; % weight function (1, C2k, sqrt) 
F = 'Gauss';  % function space 
points = 'random'; % data points (equid, Halton, Sobol, random)  

%% Prepare script 

% Test function and its integral 
CC = 20; % number of tests for Genz 
noise_level = 0; % noise added to the function values 

% initiate vectors 
N_LS = []; N_MC = []; % number of data points 
error_LS = []; error_MC = []; % errors

I = logspace(0,3,10);
%% Start some loops 
for i = 1:10 % start with d=1 for cubic PHS-RBF
            
    d = floor(I(i));
    
    %% LS RBF-LS
    % Load and store the LS-CF 
    example = matfile(['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.LS_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    error = exp_highDim( (X+1)/2, w/(2^dim), CC ); % average errors 
    % Store values
    N_LS = [N_LS; N]; % number of data points      
    error_LS = [error_LS; error]; % errors corresponding to Genz' test functions 
    
    %% (Q)MC   
    % MC weights 
    volume = 1; 
    X = (X+1)/2;
    w = ones(N,1)*volume/N;
    error = exp_highDim( X, w, CC ); % average errors 
    % Store values
    N_MC = [N_MC; N]; % number of data points      
    error_MC = [error_MC; error]; % errors corresponding to Genz' test functions 
    
end 

%% Plot figures for compparison - N vs error
figure(1) 
p = plot( N_MC,error_MC,'k^', N_LS,error_LS,'rs' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 3, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[f] - C_N[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('(Q)MC','LS','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on

