function [ f ] = interestFunction( A,P,t,x)
% this function calculates the interest function
f=A-P*(x*(1+x)^t)/((1+x)^t-1);
end

