%% generate_points
%
% Description: 
%  Function to generate the desired sample of data points
% 
% INPUT: 
%  points :     String, type of data points 
%  domain :     String, type of domain 
%  dim :        Integer, dimension 
%  omega :      weight function 
%  M :          number of initial points
%
% OUTPUT: 
%  Sample :	Structure containing 
%    dim :      the dimension 
%    volume :   the volume of the domain 
%    N :        number of data points in the domain 
%    coord :    the coordinates of the data points 
%    r :        the discrete weights needed for the LS_CF 
% 
% Author: Jan Glaubitz 
% Date: July 30, 2021 
%       

function Sample = generate_points( points, domain, dim, omega, M )

    Sample.dim = dim;

    %% Volume 
    if strcmp( domain, 'cube') 
        Sample.volume = 2^dim; % volume
    elseif strcmp( domain, 'ball') 
        Sample.volume = pi^(0.5*dim)/gamma(0.5*dim + 1); % volume 
    elseif strcmp( domain, 'combi') 
        if dim ~= 2 
           error('Desired dimension not yet implemented!') 
        end
        Sample.volume = pi + 1;
    else
        error('Desired domain not yet implemented!')
    end


    %% Let us assume a cube first
    %% equidistant point
    if strcmp( points, 'equid')  
        n = ceil( M^(1/dim) ); % number of points in every direction 
        Sample.N = n^dim; % total number of data points
        if n == 1 
            coord_aux = 0; 
        else
            coord_aux = linspace(-1, 1, n); 
        end
        % different dimensions 
        if dim == 1 
            Sample.coord = coord_aux';
        elseif dim == 2 
            [PointsX, PointsY] = meshgrid(coord_aux, coord_aux); 
            Sample.coord = [ reshape(PointsX, numel(PointsX),1)'; 
                             reshape(PointsY, numel(PointsY),1)']';
        elseif dim == 3 
            [PointsX, PointsY, PointsZ] = meshgrid(coord_aux, coord_aux, coord_aux); 
            Sample.coord = [ reshape(PointsX, numel(PointsX),1)'; 
                             reshape(PointsY, numel(PointsY),1)'; 
                             reshape(PointsZ, numel(PointsZ),1)' ]';
        else 
            error('Desired dimension not yet implemented!') 
        end   

    %% (product rule) Legendre points
    elseif strcmp( points, 'Legendre') 
        n = ceil( M^(1/dim) ); % number of points in every direction 
        Sample.N = n^dim; % total number of data points
        [x,w]=lgwt(n,-1,1);
        coord_aux = flip(x'); 
        if dim == 1
            Sample.coord = coord_aux';
        elseif dim == 2 
            [PointsX, PointsY] = meshgrid(coord_aux, coord_aux); 
            Sample.coord = [ reshape(PointsX, numel(PointsX),1)'; 
                             reshape(PointsY, numel(PointsY),1)']';
        elseif dim == 3 
            [PointsX, PointsY, PointsZ] = meshgrid(coord_aux, coord_aux, coord_aux); 
            Sample.coord = [ reshape(PointsX, numel(PointsX),1)'; 
                             reshape(PointsY, numel(PointsY),1)'; 
                             reshape(PointsZ, numel(PointsZ),1)' ]';
        else 
            error('Desired dimension not yet implemented!') 
        end

    %% random points 
    elseif strcmp( points, 'random') 
        rng('default'); rng(1,'twister'); % to make the results reproducable 
        Sample.N = M; % total number of data points
        Sample.coord = 2*rand(M,dim)-1; % generate random points

    %% Halton points 
    elseif strcmp( points, 'Halton') 
        Sample.N = M; % total number of data points
        p = haltonset(dim); % generate Halton point set 
        p = scramble(p,'RR2'); % scramble point set 
        Sample.coord = 2*net(p,M)-1; % generate the first N points

    %% Sobol points 
    elseif strcmp( points, 'Sobol') 
        Sample.N = M; % total number of data points
        p = sobolset(dim); % generate Halton point set 
        p = scramble(p,'MatousekAffineOwen'); % scramble point set 
        Sample.coord = 2*net(p,M)-1; % generate the first N points

    %% Latin hypercube sample 
    elseif strcmp( points, 'Latin') 
        Sample.N = M; % total number of data points 
        rng default % For reproducibility
        Sample.coord = 2*lhsdesign(M,dim)-1;

    %% otherwise 
    else
        error('Desired points not yet implemented!')    
    end

    %% initiate discrete weights 
    Sample.r = zeros( Sample.N, 1); % discrete weights
    Sample.r = generate_discreteWeights( Sample, omega);

    %% Cube 
    if strcmp( domain, 'cube') 
        n = 1; 
        while n <= Sample.N 
            if Sample.r(n) == 0 % weight function is zero at the points
                Sample.coord(n,:) = []; % remove point 
                Sample.r(n) = [];
                Sample.N = Sample.N - 1; % decrease number of points by 1
            else
                n = n+1; % check next point
            end 
        end

    %% Ball (only use points with radius less or equal to 1)
    elseif strcmp( domain, 'ball') 
        n = 1; 
        while n <= Sample.N 
            if norm(Sample.coord(n,:)) > 1 || Sample.r(n) == 0 % point lies outside of the ball or weight function is zero there
                Sample.coord(n,:) = []; % remove point 
                Sample.r(n) = [];
                Sample.N = Sample.N - 1; % decrease number of points by 1
            else
                n = n+1; % check next point
            end 
        end

    %% Combi (this is a combination of a ball with radius 1 with center (0,0) and a cube with radius 0.5 and center (1.5,1.5). Implemented only in 2 dimensions!)
    elseif strcmp( domain, 'combi') 
        Sample.coord = 2*Sample.coord; % stretches the data points to cover the cube with radius 2 and center (0,0) 
        n = 1; 
        while n <= Sample.N % check of points lies insider the integration domain 
            if norm(Sample.coord(n,:),2) > 1 && norm(Sample.coord(n,:)-[1.5,1.5],inf) > 0.5 % point lies outside of the domain 
                Sample.coord(n,:) = []; % remove point 
                Sample.r(n) = [];
                Sample.N = Sample.N - 1; % decrease number of points by 1
            else
                n = n+1; % check next point
            end 
        end

    %% Otherwise
    else
        error('Desired domain not yet implemented!')
    end 
    
end