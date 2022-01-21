%% generate_discreteWeights
% 
% Description: 
%  Generates the discrete weights later used to construct the LS-CF 
%
% Author: Jan Glaubitz 
% Date: July 30, 2021 
% 
% INPUT: 
%  Sample :     Sample of data points
%  omega :      weight function 
%
% OUTPUT: 
%  r         : vector of discrete weights

function r = generate_discreteWeights( Sample, omega )

    % evaluate the weight function at the data points 
    omega_evaluated = zeros(Sample.N,1);
    for n = 1:Sample.N 
       if Sample.dim == 1  
            omega_evaluated(n) = omega( Sample.coord(n,1) ); 
       elseif Sample.dim == 2  
            omega_evaluated(n) = omega( Sample.coord(n,1), Sample.coord(n,2) );
       elseif Sample.dim == 3  
            omega_evaluated(n) = omega( Sample.coord(n,1), Sample.coord(n,2) , Sample.coord(n,3) );
       elseif Sample.dim == 4  
            omega_evaluated(n) = omega( Sample.coord(n,1), Sample.coord(n,2) , Sample.coord(n,3) , Sample.coord(n,4) ); 
       elseif Sample.dim == 10  
            omega_evaluated(n) = 1;
       else 
            error('Desired dimension not yet implemented!') 
       end
    end

    % generate the discrete weights 
    r = omega_evaluated*Sample.volume/Sample.N;
    
end