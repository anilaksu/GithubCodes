% this script takes numerical integral

% Problem 1
clc;                                               % Clears the screen
clear all;

h=0.05;                                             % step size
x = 1:h:2;                                         % Calculates upto y(3)
y_e = zeros(1,length(x)); 
y_rk = zeros(1,length(x));
y_crk = zeros(1,length(x));
y_e(1) = 1;                                          % initial condition
y_rk(1) = 1;                                          % initial condition
y_crk(1) = 1;                                          % initial condition
F_xy = @(x_1,x_2) x_1+2.*x_1*x_2-1;                    % change the function as you desire

% a) Euler's method
for i=1:(length(x)-1)                              % calculation loop
    k_1 = F_xy(x(i),y_e(i));
    y_e(i+1) = y_e(i) + k_1*h;  % main equation
end
% b) Runge Kutta order 2
for i=1:(length(x)-1)                              % calculation loop
    k_1 = F_xy(x(i),y_rk(i));
    k_2 = F_xy(x(i)+0.5*h,y_rk(i)+0.5*h*k_1);
    y_rk(i+1) = y_rk(i) + (k_1+2.*k_2)*h/3;  % main equation
end
% c) Classical Runge Kutta
for i=1:(length(x)-1)                              % calculation loop
    k_1 = F_xy(x(i),y_crk(i));
    k_2 = F_xy(x(i)+0.5*h,y_crk(i)+0.5*h*k_1);
    k_3 = F_xy((x(i)+0.5*h),(y_crk(i)+0.5*h*k_2));
    k_4 = F_xy((x(i)+h),(y_crk(i)+k_3*h));
    y_crk(i+1) = y_crk(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
end

figure(1)
plot(x,y_e,'r')
hold on
plot(x,y_rk,'g')
hold on
plot(x,y_crk)
xlabel('x')
ylabel('y')
legend('Euler','2nd Order Runge Kutta', 'Classical Runge Kutta')