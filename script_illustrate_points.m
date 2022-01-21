%% script_illustrate_points  
%
% Description: 
% Script to illustrate different data points in two dimensions
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
dim = 2; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball) 
points = 'Halton'; % points (Halton, Sobol, random)

M = 8^dim; % number of points 

%% Generate and illustrate the data points 
weightFun = '1'; 
omega = generate_weightFun( weightFun, dim);
Sample = generate_points( points, domain, dim, omega, M );

figure(1)
axis equal
f1 = plot( Sample.coord(:,1), Sample.coord(:,2), 'ro', 'MarkerSize', 8); % plot points
set(f1, 'markerfacecolor', get(f1, 'color')); % use same color to fill in markers
set(gca,'FontSize',22)
