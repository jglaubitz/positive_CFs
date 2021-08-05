%% script_accuracy_Genz_compareSpaces  
%
% Description: 
%  Script to compare the accuracy of LS- and interpolatory CFs based on 
%  different function spaces for Genz' test functions
%  
% Author: Jan Glaubitz 
% Date: July 30, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
dim = 2; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball, combi) 
weightFun = '1'; % weight function (1, C2k, sqrt)
points = 'Halton'; % data points (Halton, Sobol, random) 
subsampling = 'Steinitz'; % subsampling method (NNLS, Steinitz, BP)
CC = 20; % number of tests for Genz 
noise_level = 0; % noise added to the function values 

%% Prepare some vectors 
N_LS_poly = []; % number of data points, LS-CF, polynomials 
N_LS_trig = []; % number of data points, LS-CF, trigonometric
N_LS_RBF = []; % number of data points, LS-CF, RBF
N_interpol_poly = []; % number of data points, interpolatory CF, polynomials
N_interpol_trig = []; % number of data points, interpolatory CF, trigonometric
N_interpol_RBF = []; % number of data points, interpolatory CF, RBF
error_LS_poly = []; % error, LS-CF, polynomials 
error_LS_trig = []; % error, LS-CF, trigonometric
error_LS_RBF = []; % error, LS-CF, RBF
error_interpol_poly = []; % error, interpolatory CF, polynomials
error_interpol_trig = []; % error, interpolatory CF, trigonometric
error_interpol_RBF = []; % error, interpolatory CF, RBF

% Loop over the total degree d 
for d = 1:20  
      
    [d, 20] % output to show the progress 
    
    %%% 1) Algebraic Polynomilas 
    F = 'algebraic'; 
    
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
    N_LS_poly = [N_LS_poly; N]; % number of data points      
    error_LS_poly = [error_LS_poly; error]; % errors corresponding to Genz' test functions 
    
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
    N_interpol_poly = [N_interpol_poly; N]; % number of new data points 
    error_interpol_poly = [error_interpol_poly; error]; % errors corresponding to Genz' test functions 
    
    
    %%% 2) Algebraic Polynomilas 
    F = 'trig'; 
    
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
    N_LS_trig = [N_LS_trig; N]; % number of data points      
    error_LS_trig = [error_LS_trig; error]; % errors corresponding to Genz' test functions 
    
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
    N_interpol_trig = [N_interpol_trig; N]; % number of new data points 
    error_interpol_trig = [error_interpol_trig; error]; % errors corresponding to Genz' test functions 
    
    
    %%% 3) Cubic PHS-RBF
    F = 'cubic';  
    
    %% Basis, moments, and coefficient matrix 
    omega = generate_weightFun( weightFun, dim); % set up weight function 
    [ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d, points ); % setup basis and moments  
    
    %% Load and store the LS-CF 
    example = matfile(['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.LS_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    error = Genz( (X+1)/2, w/(2^dim), CC, noise_level ); % average errors 
    % Store values
    N_LS_RBF = [N_LS_RBF; N]; % number of data points      
    error_LS_RBF = [error_LS_RBF; error]; % errors corresponding to Genz' test functions 
    
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
    N_interpol_RBF = [N_interpol_RBF; N]; % number of new data points 
    error_interpol_RBF = [error_interpol_RBF; error]; % errors corresponding to Genz' test functions 
    
end 

%% Plot figures for compparison - LS (Genz 0) 
figure(1) 
p = plot( N_LS_poly,error_LS_poly(:,1),'k^', ... 
    N_LS_trig,error_LS_trig(:,1),'rs', ... 
    N_LS_RBF,error_LS_RBF(:,1),'b+');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(3.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_0] - C_N[g_0]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('LS, poly.','LS, trig.','LS, cubic PHS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz0_compareSapces_LS_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz0_compareSapces_LS_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - interpol (Genz 0) 
figure(2) 
p = plot( N_interpol_poly,error_interpol_poly(:,1),'k^', ... 
    N_interpol_trig,error_interpol_trig(:,1),'rs', ... 
    N_interpol_RBF,error_interpol_RBF(:,1),'b+');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_0] - C_N[g_0]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('interpol., poly.','interpol., trig.','interpol., cubic PHS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz0_compareSapces_interpol_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz0_compareSapces_interpol_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - LS (Genz 1) 
figure(3) 
p = plot( N_LS_poly,error_LS_poly(:,2),'k^', ... 
    N_LS_trig,error_LS_trig(:,2),'rs', ... 
    N_LS_RBF,error_LS_RBF(:,2),'b+');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(3.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_1] - C_N[g_1]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('LS, poly.','LS, trig.','LS, cubic PHS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz1_compareSapces_LS_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz1_compareSapces_LS_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - interpol (Genz 1) 
figure(4) 
p = plot( N_interpol_poly,error_interpol_poly(:,2),'k^', ... 
    N_interpol_trig,error_interpol_trig(:,2),'rs', ... 
    N_interpol_RBF,error_interpol_RBF(:,2),'b+');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_1] - C_N[g_1]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('interpol., poly.','interpol., trig.','interpol., cubic PHS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz1_compareSapces_interpol_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz1_compareSapces_interpol_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - LS (Genz 2) 
figure(5) 
p = plot( N_LS_poly,error_LS_poly(:,3),'k^', ... 
    N_LS_trig,error_LS_trig(:,3),'rs', ... 
    N_LS_RBF,error_LS_RBF(:,3),'b+');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(3.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_2] - C_N[g_2]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('LS, poly.','LS, trig.','LS, cubic PHS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz2_compareSapces_LS_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz2_compareSapces_LS_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - interpol (Genz 2) 
figure(6) 
p = plot( N_interpol_poly,error_interpol_poly(:,3),'k^', ... 
    N_interpol_trig,error_interpol_trig(:,3),'rs', ... 
    N_interpol_RBF,error_interpol_RBF(:,3),'b+');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_2] - C_N[g_2]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('interpol., poly.','interpol., trig.','interpol., cubic PHS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz2_compareSapces_interpol_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz2_compareSapces_interpol_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - LS (Genz 3) 
figure(7) 
p = plot( N_LS_poly,error_LS_poly(:,4),'k^', ... 
    N_LS_trig,error_LS_trig(:,4),'rs', ... 
    N_LS_RBF,error_LS_RBF(:,4),'b+');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(3.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_3] - C_N[g_3]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('LS, poly.','LS, trig.','LS, cubic PHS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz3_compareSapces_LS_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz3_compareSapces_LS_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - interpol (Genz 3) 
figure(8) 
p = plot( N_interpol_poly,error_interpol_poly(:,4),'k^', ... 
    N_interpol_trig,error_interpol_trig(:,4),'rs', ... 
    N_interpol_RBF,error_interpol_RBF(:,4),'b+');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_3] - C_N[g_3]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('interpol., poly.','interpol., trig.','interpol., cubic PHS-RBF','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_Genz3_compareSapces_interpol_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_Genz3_compareSapces_interpol_dim',num2str(dim),'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 