%% script_accuracy_misc2  
%
% Description: 
%  Script to compare the accuracy of different CFs in some misc test cases 
%  
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
points = 'Halton'; % data points (Halton, Sobol, random) 

%% Fixed parameters
dim = 2; % dimension (1,2,3)
domain = 'ball'; % domain (cube, ball, combi) 
weightFun = '1'; % weight function (1, C2k, sqrt)
F = 'algebraic'; % vecor space F (algebraic, trig, RBF) 
subsampling = 'Steinitz'; % subsampling method (NNLS, Steinitz, BP)

%% Test function and reference integral 
f = @(x,y) 1./(1 + x.^2 + y.^2 ) + sin(x); % test function 
I_aux = 0.5*log(2); 
I = I_aux*( 2*pi^(dim/2)/gamma(dim/2) ); % exact integral 

%% Prepare some vectors 
N_LS = []; % number of the original data points used by the LS-CF 
N_interpol = []; % number of the subsampled data points  
N_MC = []; % number of data points used by (Q)MC 
N_Leg = []; % number of data points used by the product Legendre rule 
error_LS = []; % errors corresponding to the LS-CF 
error_interpol = []; % errors corresponding to the interpolatory CF 
error_MC = []; % errors corresponding to the (Q)MC method
error_Leg = []; % errors corresponding to the product Legendre rule

% Loop over the total degree d 
for d = 0:20 
      
    [d, 20] % output to show the progress 
    
    %% Basis, moments, and coefficient matrix 
    omega = generate_weightFun( weightFun, dim); % set up weight function 
    [ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d ); % setup basis and moments  
    
    %% Load and store the LS-CF 
    example = matfile(['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.LS_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    fun_values = f(X(:,1),X(:,2)); % function values
    error = abs( I - dot(fun_values,w) ); % average errors 
    % Store values
    N_LS = [N_LS; N]; % number of data points      
    error_LS = [error_LS; error]; % errors corresponding to Genz' test functions 
    
    %% (Q)MC   
    % MC weights 
    w = pi*omega( X(:,1), X(:,2) )/N;
    fun_values = f(X(:,1),X(:,2)); % function values 
    error = abs( I - dot(fun_values,w) ); % average errors 
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
    fun_values = f(X(:,1),X(:,2)); % function values
    error = abs( I - dot(fun_values,w) ); % average errors 
    % Store values
    N_interpol = [N_interpol; N]; % number of new data points 
    error_interpol = [error_interpol; error]; % errors corresponding to Genz' test functions 
    
    %% Load and store the product Legendre rule 
    example = matfile(['CFs/Leg_CF_dim',num2str(dim),'_',domain,'_d',num2str(d),'.mat']);
    C = example.Leg_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    fun_values = f(X(:,1),X(:,2)).*omega(X(:,1),X(:,2)); % function values 
    error = abs( I - dot(fun_values,w) ); % average errors 
    % Store values
    N_Leg = [N_Leg; N]; % number of new data points 
    error_Leg = [error_Leg; error]; % errors corresponding to Genz' test functions 
    
end 

%% Plot results 
figure(1) 
p = plot( N_MC,error_MC(:,1),'k^', ...
    N_LS,error_LS(:,1),'rs', ... 
    N_interpol,error_interpol(:,1),'b+', ... 
    N_Leg,error_Leg(:,1),'mo');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
%xlim([ 10^0, 10^3 ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[f] - C_N[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('(Q)MC','LS','interpol','Legendre','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_misc2_dim',num2str(dim),'_',domain,'_',weightFun,'_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_misc2_dim',num2str(dim),'_',domain,'_',weightFun,'_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 