function [ s ] = functionThree( x,t)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
s=0;
for i=1:t
   s=s+((-1)^(i-1))*((x^(2*i-1))/prod(2*i-1));
end


end

