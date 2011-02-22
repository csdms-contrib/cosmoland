%Modelling Landscapes Assignment 2
%Nate Bradley
%2/05/05
%returns an powerlaw distribution of random numbers


function P = powerlaw_distribution(length, gamma, x0, x1)
    if (exist('length') == 0 || ...
        exist('gamma') == 0)
        error('Not enough input arguments. Usage: powerlaw_distribution(length, gamma, lower bound, upper bound)');
    end
    
%     x0 = 10;%  10^-9; %lower bound
%     x1 = 1000000;%  1; %upper bound    
  

   
    %see http://mathworld.wolfram.com/RandomNumber.html equation 6
	P = ((x1^(gamma+1) - x0^(gamma+1))*uniform_distribution(length, 0,1)...
            + x0^(gamma+1)).^(1/(gamma+1));
	
