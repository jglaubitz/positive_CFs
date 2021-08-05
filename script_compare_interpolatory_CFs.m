%% script_compare_interpolatory_CFs  
%
% Description: 
%  Script to compare the different approaches to construct interpolatory CFs based on a given LS-CF 
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
points = 'Halton'; % data points (Halton, Sobol, Latin, random)  
CC = 20; % number of tests for Genz 
noise_level = 0; % noise added to the function values 

%% Prepare some vectors 
N_LS = []; % number of the original data points used by the LS-CF 
K = []; % degree of the function space for which th CF is exact 
N_Steinitz = []; % number of the subsampled data points produced by Steinitz 
N_BP = []; % number of the subsampled data points produced by basis pursuit
N_NNLS = []; % number of the subsampled data points produced by NNLS  
time_Steinitz = []; % time it took the Steinitz procedure 
time_BP = []; % time it took the basis pursuit approach
time_NNLS = []; % time it took the NNLS method 
error_Steinitz = []; % exactness error the Steinitz procedure 
error_BP = []; % exactness error the basis pursuit approach
error_NNLS = []; % exactness error the NNLS method 
stab_Steinitz = []; % stability measure for the Steinitz procedure 
stab_BP = []; % stability measure for the basis pursuit approach
stab_NNLS = []; % stability measure for the NNLS method 
error_Genz_LS = []; % errors corresponding to Genz' test functions 
error_Genz_Steinitz = []; % errors corresponding to Genz' test functions 
error_Genz_BP = []; % errors corresponding to Genz' test functions 
error_Genz_NNLS = []; % errors corresponding to Genz' test functions 

% Loop over the total degree d 
for d = 0:14  
      
    [d, 14] % output to show the progress 
    
    %% Basis, moments, and coefficient matrix 
    [ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d ); % setup basis and moments  
    K = [K; length(m)]; % dimension 
    
    %% Load and store the LS-CF 
    example = matfile(['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.LS_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    error_Genz = Genz( (X+1)/2, w/4, CC, noise_level ); % average errors 
    % Store values
    N_LS = [N_LS; N]; % number of data points      
    error_Genz_LS = [error_Genz_LS; error_Genz]; % errors corresponding to Genz' test functions 
    
    %% Load and store the interpolatory CF by Steinitz' method  
    example = matfile(['CFs/interpol_CF_Steinitz_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.interpol_CF_Steinitz; 
    [ N, aux] = size(C); % number of data points 
    time = C(1,dim+3); % time 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    Phi = basis(X); % coefficient matrix
    error = norm(Phi*w-m); % exactness error 
    stab = min(w); % stability measure 
    error_Genz = Genz( (X+1)/2, w/4, CC, noise_level ); % average errors 
    % Store values
    N_Steinitz = [N_Steinitz; N]; % number of new data points 
    time_Steinitz = [time_Steinitz; time]; % time 
    error_Steinitz = [error_Steinitz; error]; % exactness error 
    stab_Steinitz = [stab_Steinitz; stab]; % stability measure
    error_Genz_Steinitz = [error_Genz_Steinitz; error_Genz]; % errors corresponding to Genz' test functions 
    
    %% Load and store the interpolatory CF by basis ppursuit 
    example = matfile(['CFs/interpol_CF_BP_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.interpol_CF_BP; 
    [ N, aux] = size(C); % number of data points 
    time = C(1,dim+3); % time 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    Phi = basis(X); % coefficient matrix
    error = norm(Phi*w-m); % exactness error 
    stab = min(w); % stability measure 
    error_Genz = Genz( (X+1)/2, w/4, CC, noise_level ); % average errors 
    % Store values
    N_BP = [N_BP; N]; % number of new data points 
    time_BP = [time_BP; time]; % time 
    error_BP = [error_BP; error]; % exactness error 
    stab_BP = [stab_BP; stab]; % stability measure
    error_Genz_BP = [error_Genz_BP; error_Genz]; % errors corresponding to Genz' test functions
    
    %% Load and store the interpolatory CF by NNLS
    example = matfile(['CFs/interpol_CF_NNLS_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.interpol_CF_NNLS; 
    [ N, aux] = size(C); % number of data points 
    time = C(1,dim+3); % time 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    Phi = basis(X); % coefficient matrix
    error = norm(Phi*w-m); % exactness error 
    stab = min(w); % stability measure 
    error_Genz = Genz( (X+1)/2, w/4, CC, noise_level ); % average errors 
    % Store values
    N_NNLS = [N_NNLS; N]; % number of new data points 
    time_NNLS = [time_NNLS; time]; % time 
    error_NNLS = [error_NNLS; error]; % exactness error 
    stab_NNLS = [stab_NNLS; stab]; % stability measure
    error_Genz_NNLS = [error_Genz_NNLS; error_Genz]; % errors corresponding to Genz' test functions
    
end 

%% Plot figures for compparison - K vs N 
figure(1) 
p = plot( K,N_Steinitz,'rs', K,N_BP,'b+', K,N_NNLS,'k^' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ max([NN_Leg(1);NN_LS(1)]), min([NN_Leg(end);NN_LS(end)]) ]) 
%ylim([ , ])
xlabel('$K$','Interpreter','latex') 
ylabel('$N$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('Steinitz','basis pursuit','NNLS','Location','southeast'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/compare_interpolCFs_KvsN_dim',num2str(dim),'_',points,'.fig'] );
%savefig(str);
str = sprintf( ['plots/compare_interpolCFs_KvsN_dim',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - N_LS vs time 
% normalize the times
time_max = max([time_Steinitz;time_BP;time_NNLS]); 
time_Steinitz = time_Steinitz/time_max; 
time_BP = time_BP/time_max; 
time_NNLS = time_NNLS/time_max; 
% plot figure 
figure(2) 
p = plot( N_LS,time_Steinitz,'rs', N_LS,time_BP,'b+', N_LS,time_NNLS,'k^' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ max([NN_Leg(1);NN_LS(1)]), min([NN_Leg(end);NN_LS(end)]) ]) 
%ylim([ , ])
xlabel('$N$','Interpreter','latex') 
ylabel('time','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('Steinitz','basis pursuit','NNLS','Location','northwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/compare_interpolCFs_NvsTime_dim',num2str(dim),'_',points,'.fig'] );
%savefig(str);
str = sprintf( ['plots/compare_interpolCFs_NvsTime_dim',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - N vs error 
figure(3) 
p = plot( N_LS,error_Steinitz,'rs', N_LS,error_BP,'b+', N_LS,error_NNLS,'k^' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ max([NN_Leg(1);NN_LS(1)]), min([NN_Leg(end);NN_LS(end)]) ]) 
%ylim([ , ])
xlabel('$N$','Interpreter','latex') 
ylabel('$\| \Phi \mathbf{w} - \mathbf{m} \|_2$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('Steinitz','basis pursuit','NNLS','Location','northwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/compare_interpolCFs_NvsError_dim',num2str(dim),'_',points,'.fig'] );
%savefig(str);
str = sprintf( ['plots/compare_interpolCFs_NvsError_dim',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - N vs stability 
figure(4) 
p = plot( N_LS,stab_Steinitz,'rs', N_LS,stab_BP,'b+', N_LS,stab_NNLS,'k^' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ max([NN_Leg(1);NN_LS(1)]), min([NN_Leg(end);NN_LS(end)]) ]) 
%ylim([ , ])
xlabel('$N$','Interpreter','latex') 
ylabel('$\min \mathbf{w}$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('Steinitz','basis pursuit','NNLS','Location','northeast'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/compare_interpolCFs_NvsStab_dim',num2str(dim),'_',points,'.fig'] );
%savefig(str);
str = sprintf( ['plots/compare_interpolCFs_NvsStab_dim',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - N vs error (Genz 0) 
figure(5) 
p = plot( N_Steinitz,error_Genz_Steinitz(:,1),'rs', N_BP,error_Genz_BP(:,1),'b+', ... 
    N_NNLS,error_Genz_NNLS(:,1),'k^', N_LS,error_Genz_LS(:,1),'mo' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^2 ]) 
%ylim([ , ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_1] - C_N[g_1]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('Steinitz','basis pursuit','NNLS','LS','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/compare_interpolCFs_Genz0_dim',num2str(dim),'_',points,'.fig'] );
%savefig(str);
str = sprintf( ['plots/compare_interpolCFs_Genz0_dim',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - N vs error (Genz 1) 
figure(6) 
p = plot( N_Steinitz,error_Genz_Steinitz(:,2),'rs', N_BP,error_Genz_BP(:,2),'b+', ... 
    N_NNLS,error_Genz_NNLS(:,2),'k^', N_LS,error_Genz_LS(:,2),'mo' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^2 ]) 
%ylim([ , ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_1] - C_N[g_1]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('Steinitz','basis pursuit','NNLS','LS','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/compare_interpolCFs_Genz1_dim',num2str(dim),'_',points,'.fig'] );
%savefig(str);
str = sprintf( ['plots/compare_interpolCFs_Genz1_dim',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - N vs error (Genz 2) 
figure(7) 
p = plot( N_Steinitz,error_Genz_Steinitz(:,3),'rs', N_BP,error_Genz_BP(:,3),'b+', ... 
    N_NNLS,error_Genz_NNLS(:,3),'k^', N_LS,error_Genz_LS(:,3),'mo' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ 10^0, 10^2 ]) 
%ylim([ , ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_2] - C_N[g_2]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('Steinitz','basis pursuit','NNLS','LS','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/compare_interpolCFs_Genz2_dim',num2str(dim),'_',points,'.fig'] );
%savefig(str);
str = sprintf( ['plots/compare_interpolCFs_Genz2_dim',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc') 

%% Plot figures for compparison - N vs error (Genz 3) 
figure(8) 
p = plot( N_Steinitz,error_Genz_Steinitz(:,4),'rs', N_BP,error_Genz_BP(:,4),'b+', ... 
    N_NNLS,error_Genz_NNLS(:,4),'k^', N_LS,error_Genz_LS(:,1),'mo' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^2 ]) 
%ylim([ , ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[g_3] - C_N[g_3]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('Steinitz','basis pursuit','NNLS','LS','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/compare_interpolCFs_Genz3_dim',num2str(dim),'_',points,'.fig'] );
%savefig(str);
str = sprintf( ['plots/compare_interpolCFs_Genz3_dim',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc') 