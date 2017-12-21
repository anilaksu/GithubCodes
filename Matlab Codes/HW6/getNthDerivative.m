function [ f_n ] = getNthDerivative( n,x_0 )
% this function calculates the n-th derivative of the function
f_n=factorial(n)*(1-x_0)^(-n-1);
end

