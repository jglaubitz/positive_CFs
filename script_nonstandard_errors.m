%% script_nonstandard_errors 
%
% Description: 
% Script to compare the interpolatory and LS-CFs for a nonstandard domain with QMC  
%
% Author: Jan Glaubitz 
% Date: Jan 13, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
points = 'random'; % data points (equid, Halton, Sobol, Latin, random)  
d_max = 14; % max degree of exactness 

% Fixed parameters
dim = 2; % dimension 
F = 'algebraic'; % function space 
weightFun = '1'; % weight function 
domain = 'nonstandard'; 


%% Prepare script 

% Test function and its integral 
fun_test = @(x,y) exp(x.^2 + y.^2); 
fun_r = @(theta) 1 - (1/3)*sin(2*theta).^2; % parameterization of boundary  
fun_int = @(theta) exp(fun_r(theta).^2); 
int_ref = 0.5*integral(fun_int,0,2*pi) - pi;

% initiate vectors 
N_LS = []; N_interpol = []; N_MC = []; % number of data points 
error_LS = []; error_interpol = []; error_MC = []; % errors 

% Loop over the total degree d 
for d = 0:d_max % start with d=1 for cubic PHS-RBF
            
    %% LS-CF 
    example = matfile(['CFs/LS_CF_',domain,'_d',num2str(d),'_',points,'.mat']);
    C = example.LS_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    fun_values = fun_test(X(:,1),X(:,2)); % function values 
    error = abs( int_ref - dot(fun_values,w) ); % average errors 
    % Store values
    N_LS = [N_LS; N]; % number of data points      
    error_LS = [error_LS; error]; % errors corresponding to Genz' test functions 
    
    %% (Q)MC   
    % MC weights 
    volume = 17*pi/24; 
    w = ones(N,1)*volume/N;
    fun_values = fun_test(X(:,1),X(:,2)); % function values 
    error = abs( int_ref - dot(fun_values,w) ); % average errors 
    % Store values
    N_MC = [N_MC; N]; % number of data points      
    error_MC = [error_MC; error]; % errors corresponding to Genz' test functions 
    
     %% interpolatory CF 
    example = matfile(['CFs/interpol_CF_',domain,'_d',num2str(d),'_',points,'.mat']);
    C = example.interpol_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    fun_values = fun_test(X(:,1),X(:,2)); % function values 
    error = abs( int_ref - dot(fun_values,w) ); % average errors 
    % Store values
    N_interpol = [N_interpol; N]; % number of data points      
    error_interpol = [error_interpol; error]; % errors corresponding to Genz' test functions 
    
end

%% Plot results 
figure(1) 
p = plot( N_MC,error_MC(:,1),'k^', ...
    N_LS,error_LS(:,1),'rs', ... 
    N_interpol,error_interpol(:,1),'b+');
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ 10^0, 10^3 ]) 
%ylim([ 10^(-6), 10^1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[f] - C_N[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('(Q)MC','LS','interpol','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on