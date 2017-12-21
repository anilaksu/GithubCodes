% this script is developed to solve HW5

%Problem 1
x=linspace(2,6,5);

f=4.*x.^3-2;

figure(1)

plot(x,f)
xlabel('x', 'FontSize', 30)
ylabel('f(x)', 'FontSize', 30)

% let's calculate the trapezoids under the curve

for i=1:(length(x)-1)
   A(i)=0.5*(f(i)+f(i+1)); 
end

% Problem 2
% the tolerance
tol=1E4;
% the first indice
n=1;
% first element of Woodall
W(1)=1;
 while 1
   n=n+1;
   W_aux=n*2^n-1;
   if W_aux > tol
       break
   end
   W(n)=W_aux
 end

% the sum of elements
sum_W=0.;
for i=1:length(W)
   sum_W=sum_W+W(i); 
end
% the average of Woodall Numbers
WA=sum_W/length(W);

% Problem 3
t=linspace(0,4,21);
% the initial velocity
V_0=20.;
% the initial height 
y_0=20.;
% the initial angle
theta=47;
% the gravity constant
g=9.81;
% the y locations
y=y_0+V_0*sind(theta)*t-0.5*g*t.^2;
% the x locations
x=V_0*cosd(theta)*t;

figure(2)

plot(x,y)
xlabel('x', 'FontSize', 30)
ylabel('y', 'FontSize', 30)

% Problem 4
f(1)=1;
f(2)=1;
for i=3:10
    f(i)=f(i-1)+f(i-2);
end

% Problem 5

% the order of the derivative
n=5;
% the location
x=2;
% the n-th derivative of function at x
f=((-1)^n)*factorial(n+1)/(x^(n+2));

% Problem 6
% the base axis 
b=linspace(400,600,11);

% the inertia
A_1=40.*400.;
A_2=40.*320.;
A_3=40.*b;
% the neutral axis
for i=1:length(b)
    yd(i)=(A_3(i)*20+A_2*200+A_1*380)/(A_1+A_2+A_3(i));
    I_1(i)=A_1*(380-yd(i))^2+400*(40^3)/12;
    I_2(i)=A_2*(200-yd(i))^2+40*(320^3)/12;
    I_3(i)=A_3(i)*(20-yd(i))^2+b(i)*(40^3)/12;
    I_tot(i)=I_1(i)+I_2(i)+I_3(i);
end

% Problem 7 
tol=1E-5;
% the initial guess
a(1)=1;
% the initial error
err_7(1)=1.;
% indice
i=1;
while err_7(i) > tol
   i=i+1;
   a(i)=0.5*a(i-1)+1./a(i-1);
   err_7(i)=abs((a(i)-a(i-1))/a(i));
end

% Problem 8
P=40000;
A=12000;
t=5;
% the tolerance 
tol=1E-2;
% initial err
err_8=1.;
% the right point
x_r=0.2;
% the left point
x_l=0.05;
while err_8 > tol
    % the function at x_l
    f_r = interestFunction( A,P,t,x_r);
    % the function at x_r
    f_l = interestFunction( A,P,t,x_l);
    % the mid point
    x_m=0.5*(x_r+x_l);
    % the function at x_m
    f_m= interestFunction( A,P,t,x_m);
    % the error 
    err_8 = abs(f_m);
    
    if(f_r*f_m > 0)
        x_r=x_m;
    else
        x_l=x_m;
    end
    
end
