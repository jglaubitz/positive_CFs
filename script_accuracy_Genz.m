%% script_accuracy_Genz  
%
% Description: 
%  Script to compare the accuracy of different CFs for Genz' test functions
%  
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
dim = 2; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball, combi) 
weightFun = '1'; % weight function (1, C2k, sqrt)
F = 'algebraic'; % vecor space F (algebraic, trig, RBF) 
points = 'Halton'; % data points (Halton, Sobol, random) 
subsampling = 'Steinitz'; % subsampling method (NNLS, Steinitz, BP)
CC = 20; % number of tests for Genz 
noise_level = 0; % noise added to the function values 

%% Prepare some vectors 
N_LS = []; % number of the original data points used by the LS-CF 
N_interpol = []; % number of the subsampled data points  
N_MC = []; % number of data points used by (Q)MC 
N_Leg = []; % number of data points used by the product Legendre rule 
error_LS = []; % errors corresponding to the LS-CF 
error_interpol = []; % errors corresponding to the interpolatory CF 
error_MC = []; % errors corresponding to the (Q)MC method
error_Leg = []; % errors corresponding to the product Legendre rule

% max total degree
if dim == 2 
    d_max = 20; 
else 
    d_max = 12; 
end

% Loop over the total degree d 
for d = 0:d_max
      
    [d,d_max] % output to show progress 
    
    %% Basis, moments, and coefficient matrix 
    omega = generate_weightFun( weightFun, dim); % set up weight function 
    [ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d ); % setup basis and moments  
    
    %% Load and store the LS-CF 
    example = matfile(['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.LS_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    error = Genz( (X+1)/2, w/(2^dim), CC, noise_level ); % average errors 
    % Store values
    N_LS = [N_LS; N]; % number of data points      
    error_LS = [error_LS; error]; % errors corresponding to Genz' test functions 
    
    %% (Q)MC   
    % MC weights 
    if dim == 1  
      	w = omega( X(:,1) )/N;
  	elseif dim == 2  
       	w = omega( X(:,1), X(:,2) )/N;
   	elseif dim == 3  
       	w = omega( X(:,1), X(:,2), X(:,3) )/N;
    else 
       	error('Desired dimension not yet implemented!') 
    end 
    error = Genz( (X+1)/2, w, CC, noise_level ); % average errors 
    % Store values
    N_MC = [N_MC; N]; % number of data points      
    error_MC = [error_MC; error]; % errors corresponding to Genz' test functions 
    
    %% Load and store the interpolatory CF 
    example = matfile(['CFs/interpol_CF_',subsampling,'_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);  
    if strcmp( subsampling, 'Steinitz')
        C = example.interpol_CF_Steinitz; 
    elseif strcmp( subsampling, 'NNLS')
        C = example.interpol_CF_NNLS;
    elseif strcmp( subsampling, 'BP')
        C = example.interpol_CF_BP;
    else 
        error('Unknown subsampling method!')
    end
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    error = Genz( (X+1)/2, w/(2^dim), CC, noise_level ); % average errors 
    % Store values
    N_interpol = [N_interpol; N]; % number of new data points 
    error_interpol = [error_interpol; error]; % errors corresponding to Genz' test functions 
    
    %% Load and store the product Legendre rule 
    example = matfile(['CFs/Leg_CF_dim',num2str(dim),'_',domain,'_d',num2str(d),'.mat']);
    C = example.Leg_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    error = Genz( (X+1)/2, w/(2^dim), CC, noise_level ); % average errors 
    % Store values
    N_Leg = [N_Leg; N]; % number of new data points 
    error_Leg = [error_Leg; error]; % errors corresponding to Genz' test functions 
    
end 

%% Plot figures for compparison - N vs error (Genz 0) 
figure(1) 
p = plot( N_MC,error_MC(:,1),'k^', ... 
    N_LS,error_LS(:,1),'rs', ... 
    N_interpol,error_interpol(:,1),'b+', ... 
    N_Leg,error_Leg(:,1),'mo');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(3.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_0] - C_N[g_0]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('(Q)MC','LS','interpol','Legendre','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz0_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz0_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - N vs error (Genz 1) 
figure(2) 
p = plot( N_MC,error_MC(:,2),'k^', ...
    N_LS,error_LS(:,2),'rs', ... 
    N_interpol,error_interpol(:,2),'b+', ... 
    N_Leg,error_Leg(:,2),'mo');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(3.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_1] - C_N[g_1]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('(Q)MC','LS','interpol','Legendre','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz1_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz1_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - N vs error (Genz 2) 
figure(3) 
p = plot( N_MC,error_MC(:,3),'k^', ...
    N_LS,error_LS(:,3),'rs', ... 
    N_interpol,error_interpol(:,3),'b+', ... 
    N_Leg,error_Leg(:,3),'mo');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(3.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_2] - C_N[g_2]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('(Q)MC','LS','interpol','Legendre','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz2_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz2_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - N vs error (Genz 3) 
figure(4) 
p = plot( N_MC,error_MC(:,4),'k^', ...
    N_LS,error_LS(:,4),'rs', ... 
    N_interpol,error_interpol(:,4),'b+', ... 
    N_Leg,error_Leg(:,4),'mo');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(3.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_3] - C_N[g_3]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('(Q)MC','LS','interpol','Legendre','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz3_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz3_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 