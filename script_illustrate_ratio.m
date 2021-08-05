%% script_illustrate_ratio 
%
% Description: 
%  Script to illustrate the relation between the number of data points N
%  and the dimension K of the function space that is treated exactly by the
%  LS-CF 
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
points = 'Halton'; % data points (Halton, Sobol, random)  

%% Prepare some vectors 
K = []; % degree of the function space for which th CF is exact 
K_Leg = []; % degree of the function space for which th CF is exact
N_LS = []; % number of the original data points used by the LS-CF 
N_interpol = []; % number of the subsampled data points produced by NNLS  
N_Leg = []; % number of data points used by the product Legendre rule 

% Loop over the total degree d 
for d = 0:10  
      
    %% Basis, moments, and coefficient matrix 
    [ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d, points ); % setup basis and moments  
    K = [K; length(m)]; % dimension 
    
    %% Load and store the LS-CF 
    example = matfile(['CFs/LS_CF_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.LS_CF; 
    [ N, aux] = size(C); % number of data points 
    % Store values
    N_LS = [N_LS; N]; % number of data points      
       
    %% Load and store the interpolatory CF by NNLS
    example = matfile(['CFs/interpol_CF_NNLS_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_d',num2str(d),'_',points,'.mat']);
    C = example.interpol_CF_NNLS; 
    [ N, aux] = size(C); % number of data points 
    % Store values
    N_interpol = [N_interpol; N]; % number of new data points 
    
    %% Load and store the product Legendre rule 
    example = matfile(['CFs/Leg_CF_dim',num2str(dim),'_',domain,'_d',num2str(d),'.mat']);
    C = example.Leg_CF; 
    [ N, aux] = size(C); % number of data points 
    % Store values
    N_Leg = [N_Leg; N]; % number of new data points 
    K_Leg = [K_Leg; C(1,dim+3)];
    
end 

%% Plot figures for compparison - K vs N 
figure(1) 
p = plot( K,N_LS,'rs', K,N_interpol,'b+', K_Leg,N_Leg,'k^' );
set(p, 'LineWidth',2)
set(p, 'markersize',10)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ 1, 10^2]) 
%ylim([ , ])
xlabel('$K$','Interpreter','latex') 
ylabel('$N$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
lgnd = legend('LS','interpolatory','Legendre','Location','northwest'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on
str = sprintf( ['plots/ratio_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_',points,'.fig'] );
%savefig(str);
str = sprintf( ['plots/ratio_dim',num2str(dim),'_',domain,'_',weightFun,'_',F,'_',points,'.eps'] );
%saveas(gcf,str,'epsc') 

% LS fit for the parameters s and C - LS 
f_LS = fittype('C*x.^s'); % set up model for parameters 
xdata = N_LS; 
ydata = K; 
fit( ydata, xdata, f_LS, 'Lower', [0,0] ) 