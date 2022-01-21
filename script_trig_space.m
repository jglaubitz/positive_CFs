%% script_trig_space  
%
% Description: 
%  Script to compare the accuracy of LS- and interpolatory CFs for
%  polynomial and trigonometric function spaces 
%  
% Author: Jan Glaubitz 
% Date: Jan 11, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters 
domain = 'cube'; % domain (cube, ball, combi) 
weightFun = '1'; % weight function (1, C2k, sqrt)
points = 'random'; % data points (equid, Halton, Sobol, random) 
subsampling = 'Steinitz'; % subsampling method (NNLS, Steinitz, BP) 

%% Prepare script 

% Test function and its integral 
dim = 1; % dimension 
fun_test = @(x) exp(sin(pi*x)).*cos(pi*x); 
int_ref = integral(fun_test,-1,1);

% initiate vectors 
N_LS_poly = []; % number of data points, LS-CF, polynomials 
N_LS_trig = []; % number of data points, LS-CF, trigonometric
N_trap = []; % number of data points, trapezoidal rule
error_LS_poly = []; % error, LS-CF, polynomials 
error_LS_trig = []; % error, LS-CF, trigonometric
error_trap = []; % error, trapezoidal rule

% Loop over the total degree d 
for d = 1:30  
      
    [d, 30] % output to show the progress 
    
    %% 1) Algebraic Polynomilas 
    F = 'algebraic';   
    
    % Load and store the LS-CF 
    example = matfile(['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.LS_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    error = abs( int_ref - dot( w, fun_test(X) ) ); % error
    % Store values
    N_LS_poly = [N_LS_poly; N]; % number of data points      
    error_LS_poly = [error_LS_poly; error]; % errors corresponding to Genz' test functions 
    
    %% 2) Trigonometric functions 
    F = 'trig'; 
  
    % Load and store the LS-CF 
    example = matfile(['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.LS_CF; 
    [ N, aux] = size(C); % number of data points 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    error = abs( int_ref - dot( w, fun_test(X) ) ); % error 
    % Store values
    N_LS_trig = [N_LS_trig; N]; % number of data points      
    error_LS_trig = [error_LS_trig; error]; % errors corresponding to Genz' test functions 
    
    %% 3) Trapezoidal rule 
    X = linspace(-1, 1, N)'; 
    w = 2/(N-1)*ones(N,1); 
    w(1) = 0.5*w(1); w(end) = 0.5*w(end); 
    error = abs( int_ref - dot( w, fun_test(X) ) ); % error 
    N_trap = [N_trap; N]; % number of data points  
    error_trap = [error_trap; error]; % errors corresponding to Genz' test functions 
    
end 

%% Plot figures for compparison - LS 
figure(1) 
p = plot( N_LS_poly,error_LS_poly,'k^', ... 
    N_LS_trig,error_LS_trig,'rs' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ 3, 10^(2.5) ]) 
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[f] - C_N[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('LS poly','LS trig','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_trig_space_LS_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_trig_space_LS_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 

% include trapezoidal rule 
figure(2) 
p = plot( N_LS_poly, error_LS_poly, 'k^', ... 
    N_LS_trig, error_LS_trig, 'rs', ... 
    N_trap, error_trap, 'bo' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ 3, 10^(2.5) ])
%ylim([ 10^(-15), 1])
xlabel('$N$','Interpreter','latex') 
ylabel('$|I[f] - C_N[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('LS poly','LS trig','trapezoidal','Location','southwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/accuracy_trig_space_LS_',points,'_',subsampling,'.fig'] );
%savefig(str);
str = sprintf( ['plots/accuracy_trig_space_LS_',points,'_',subsampling,'.eps'] );
%saveas(gcf,str,'epsc') 