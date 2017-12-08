function [ f ] = FunctionFor( t )
% this function calculates the values of f(t)=t*cos(t)-t

for i=1:length(t)
    f(i)=t(i)*cos(t(i))-t(i);
end

end

