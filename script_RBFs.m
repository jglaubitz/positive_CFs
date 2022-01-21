%% script_RBFs
%
% Description: 
%  Script to compare the accuracy and stbility of interpolatory RBF-CFs and
%  LS RBF-CFs 
%  
% Author: Jan Glaubitz 
% Date: Jan 12, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters 
domain = 'cube'; % domain (cube, ball, combi) 
weightFun = '1'; % weight function (1, C2k, sqrt) 
F = 'Gauss';  % function space 
points = 'random'; % data points (equid, Halton, Sobol, random)  

%% Prepare script 

% Test function and its integral 
dim = 2; % dimension 
CC = 20; % number of tests for Genz 
noise_level = 2; % noise added to the function values 

% initiate vectors 
N_RBF = []; N_LS = []; % number of data points 
error_RBF = []; error_LS = []; % errors
stab_RBF = []; stab_LS = []; % stability measure
w_min_RBF = []; w_min_LS = []; % minimal weight

% Loop over the total degree d 
for k = 1:60  
    
    %% 1) RBF-CF    
    % Load and store the LS-CF 
    example = matfile(['CFs/RBF_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(k),'_',points,'.mat']);
    C = example.RBF_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    error = Genz( (X+1)/2, w/(2^dim), CC, noise_level ); % average errors 
    stab = sum(abs(w)); % stability measure 
    w_min = min(w); 
    % Store values
    N_RBF = [N_RBF; N]; % number of data points      
    error_RBF = [error_RBF; error]; % errors corresponding to Genz' test functions 
    stab_RBF = [stab_RBF; stab]; % stability measure 
    w_min_RBF = [w_min_RBF; w_min];
    
    %% 2) LS RBF-LS
    % Load and store the LS-CF 
    example = matfile(['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(k),'_',points,'.mat']);
    C = example.LS_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    error = Genz( (X+1)/2, w/(2^dim), CC, noise_level ); % average errors 
    stab = sum(abs(w)); % stability measure 
    w_min = min(w); 
    % Store values
    N_LS = [N_LS; N]; % number of data points      
    error_LS = [error_LS; error]; % errors corresponding to Genz' test functions 
    stab_LS = [stab_LS; stab]; % stability measure 
    w_min_LS = [w_min_LS; w_min];
    
end 

[w_min_RBF, w_min_LS]

%% Plot figures for compparison - N vs error (Genz 0) 
figure(1) 
p = plot( N_RBF,error_RBF(:,1),'k^', N_LS,error_LS(:,1),'rs' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 3, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[f] - C_N[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('RBF','LS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on

%% Plot figures for compparison - N vs error (Genz 1) 
figure(2) 
p = plot( N_RBF,error_RBF(:,2),'k^', N_LS,error_LS(:,2),'rs' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 3, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[f] - C_N[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('RBF','LS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on

%% Plot figures for compparison - N vs error (Genz 2) 
figure(3) 
p = plot( N_RBF,error_RBF(:,3),'k^', N_LS,error_LS(:,3),'rs' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 3, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[f] - C_N[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('RBF','LS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on

%% Plot figures for compparison - N vs error (Genz 3) 
figure(4) 
p = plot( N_RBF,error_RBF(:,4),'k^', N_LS,error_LS(:,4),'rs' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 3, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[f] - C_N[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('RBF','LS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on

%% Plot - N vs minimal weight 
[B,I1] = sort(N_RBF); [B,I2] = sort(N_LS); 
figure(5) 
p = plot( N_RBF(I1),w_min_RBF(I1),'k--', N_LS(I2),w_min_LS(I2),'r:' );
set(p, 'LineWidth',3)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 3, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$\mathrm{min}(\mathbf{w})$','Interpreter','latex')
set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
lgnd = legend('RBF','LS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on