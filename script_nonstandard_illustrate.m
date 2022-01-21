%% script_nonstandard_illustrate   
%
% Description: 
% Script to illustrate the nonstandard domain  
%
% Author: Jan Glaubitz 
% Date: Jan 13, 2021 

%% Setting up the script 
clc, clear 

%% Free parameters
points = 'random'; % data points (equid, Halton, Sobol, Latin, random)  
M = 200; % degree of exactness 

% Fixed parameters
dim = 2; % dimension 
F = 'algebraic'; % function space 
weightFun = '1'; % weight function 
domain = 'nonstandard'; 

%% Generate data points 
omega = generate_weightFun( weightFun, dim);
Sample = generate_points( points, domain, dim, omega, M );
X = Sample.coord; 

%% Compute the boundary of the domain 
nr_points = 200; 
theta = linspace(0,2*pi,nr_points); % angles 
fun_r = @(theta) 1 - (1/3)*sin(2*theta).^2; % parameterization of boundary  
r = fun_r(theta); % corresponding radii 
[x,y] = pol2cart(theta,r); % boundary in Cartesian coordinates 

%% Illustrate 
figure(1) 
hold on 
f1 = plot(x,y,'b'); 
set(f1, 'linewidth', 4); % use same color to fill in markers
f2 = plot( X(:,1), X(:,2), 'ro', 'MarkerSize', 8); % plot points
set(f2, 'markerfacecolor', get(f2, 'color')); % use same color to fill in markers
set(gca,'FontSize',22) 
xlim([-1,1]); ylim([-1,1]);
hold off