function [ f_d,g_d ] = calculateDerivative( x_d )
% this function computes the derivative of f function
syms f(x) g(x)
f(x) = 3.*x^3+2.*x^2-x+4;
g(x) = diff(f,x);

for i=1:length(x_d)
    f_d(i)=f(x_d(i));
    g_d(i)=g(x_d(i));
end

end

