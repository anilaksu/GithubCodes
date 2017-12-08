function [ g ] = FunctionWhile( t )
%this function calculates the values of g(t)=exp(t)-t

% let's give the initial value of the indice
i=1;
while i<=length(t)
    g(i)=exp(t(i))-t(i);
    i=i+1;
end

end

