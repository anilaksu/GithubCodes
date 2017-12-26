function [ y ] = getStateSpaceEuler( x,y )
%this function calculates the state space form of the right hand-side of 
% ODE given in HW 12 by Euler's method

h=abs(x(2)-x(1));                                % the step size
F_2 = @(x_1,x_2,x_3) x_1+x_2+x_3;                % change the function as you desire
F_1 = @(x_3) x_3;                                % change the function as you desire

% a) Euler's method
for i=1:(length(x)-1)                              % calculation loop
    k_1(1,1) = F_1(y(2,i));
    k_1(2,1) = F_2(x(i),y(1,i),y(2,i));    
    y(:,i+1) = y(:,i)+k_1*h;  % main equation
end

end

