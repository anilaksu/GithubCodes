% this script takes numerical time integral

% Problem 1
clc;                                               % Clears the screen
clear all;

h=0.1;                                             % step size
x = 0:h:2;                                         % Calculates upto y(2)

y_e = zeros(2,length(x)); 
y_m = zeros(2,length(x));
y_crk = zeros(2,length(x));
y_abm = zeros(2,length(x));

y_e(1,1) = 2;                                          % initial condition
y_m(1,1) = 2;                                          % initial condition
y_crk(1,1) = 2;                                        % initial condition
y_abm(1,1) = 2;                                        % initial condition

F_2 = @(x_1,x_2,x_3) x_1+x_2+x_3;                    % change the function as you desire
F_1 = @(x_3) x_3;                                % change the function as you desire

% a) Euler's method
for i=1:(length(x)-1)                              % calculation loop
    k_1(1,1) = F_1(y_e(2,i));
    k_1(2,1) = F_2(x(i),y_e(1,i),y_e(2,i));    
    y_e(:,i+1) = y_e(:,i)+k_1*h;  % main equation
end
% b) Mid Point Method
for i=1:(length(x)-1)                              % calculation loop
    k_1(1,1) = F_1(y_e(2,i));
    k_1(2,1) = F_2(x(i),y_m(1,i),y_m(2,i)); 
    k_2(1,1) = F_1(y_m(2,i)+0.5*h*k_1(2,1));
    k_2(2,1) = F_2(x(i)+0.5*h,y_m(1,i)+0.5*h*k_1(1,1),y_m(2,i)+0.5*h*k_1(2,1)); 
    y_m(:,i+1) = y_m(:,i) + (k_1+2.*k_2)*h/3;  % main equation
end
 % c) Classical Runge Kutta
 for i=1:(length(x)-1)                              % calculation loop
    k_1(1,1) = F_1(y_crk(2,i));
    k_1(2,1) = F_2(x(i),y_crk(1,i),y_crk(2,i)); 
    k_2(1,1) = F_1(y_m(2,i)+0.5*h*k_1(2,1));
    k_2(2,1) = F_2(x(i)+0.5*h,y_crk(1,i)+0.5*h*k_1(1,1),y_crk(2,i)+0.5*h*k_1(2,1)); 
    k_3(1,1) = F_1(y_m(2,i)+0.5*h*k_2(2,1));
    k_3(2,1) = F_2(x(i)+0.5*h,y_crk(1,i)+0.5*h*k_2(1,1),y_crk(2,i)+0.5*h*k_2(2,1)); 
    k_4(1,1) = F_1(y_m(2,i)+h*k_3(2,1));
    k_4(2,1) = F_2(x(i)+h,y_crk(1,i)+h*k_3(1,1),y_crk(2,i)+h*k_3(2,1)); 
    y_crk(:,i+1) = y_crk(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
 end


% d) Adams-Bashforth-Moulton Method
% in this method first 3 steps  are integrated by other methods
[ y_abm ] = getStateSpaceEuler( x(1:3),y_abm );

for i=4:(length(x))                              % calculation loop
    k_1(1,1) = F_1(y_abm(2,i-2));
    k_1(2,1) = F_2(x(i-2),y_abm(1,i-2),y_abm(2,i-2)); 
    k_2(1,1) = F_1(y_abm(2,i-1));
    k_2(2,1) = F_2(x(i-1),y_abm(1,i-1),y_abm(2,i-1)); 
    P_1(:,1) = y_abm(:,i-1) + (1/2)*(3*k_2-k_1)*h;  % predictor step
    k_3(1,1) = F_1(P_1(2,1));
    k_3(2,1) = F_2(x(i),P_1(1,1),P_1(2,1)); 
    y_abm(:,i) = y_abm(:,i-1) + (1/2)*(k_3+k_2)*h;  % corrector step
end

figure(1)
plot(x,y_e(1,:),'r')
hold on
plot(x,y_m(1,:),'g')
hold on
plot(x,y_crk(1,:))
hold on
plot(x,y_abm(1,:),'m')
xlabel('x')
ylabel('y')
legend('Euler','Mid-point', 'Classical Runge Kutta','Adams-Bashforth-Moulton')