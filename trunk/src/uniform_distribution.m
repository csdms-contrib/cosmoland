%Modelling Landscapes Assignment 2
%Nate Bradley
%2/05/05
%return a uniform distribution of random numbers

function U = uniform_distribution(length, min, max)
    if (exist('length') == 0 || ...
        exist('min') == 0 || ...
        exist('max') == 0)
        error('Not enough input arguments. Usage: uniform_distribution(length, min, max)');
    end
    U = (rand(length, 1) * (max-min)) + min;


