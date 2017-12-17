% this script is developed to test our while skills

%Problem 1
a = [1 4 9 16 25 36];

% summation of a array
sum_a=0.;

for i=1:length(a)
   sum_a=sum_a+a(i); 
end

%Problem 3
x=1.2;
y=2.3;
z=3.9;

[ r ] = getSqrtProd( x,y,z );

%Problem 4
% let's define b vector by linspace command
b=linspace(0,12,7);
% let's compute c vector
c=b.^3;

%Problem 5
for i=1:(length(c)-1)
   d(i)=0.5*(c(i)+c(i+1)); 
end

%Problem 6
for i=1:(length(c)-2)
   e(i)=(c(i)+c(i+1)+c(i+2))/3.; 
end

%Problem 7 
E=200E9;
I=1E-4;
L=1.2;
w=1000;
% let's define x vector
x=linspace(0,1.2,7);

for i=1:length(x)
    del(i)=(-w*x(i)^2.)*(6.*L^2. - 4.*L*x(i)+x(i)^2.)/(24.*E*I);
end

figure(1)

plot(x,del)
xlabel('x', 'FontSize', 30)
ylabel('\delta', 'FontSize', 30)

% Problem 8
% syms f(x) g(x)
% f(x) = 3.*x^3+2.*x^2-x+4;
% g(x) = diff(f,x);

x=linspace(4,12,9);
[ f_d,g_d ] = calculateDerivative( x )

figure(2)
plot(x,f_d)
xlabel('x', 'FontSize', 30)
ylabel('f(x)', 'FontSize', 30)

% Problem 9
% first element of the series
pi_1(1)=4.
pi_2(1)=2.
for i=2:15
    pi_1(i)=pi_1(i-1)+4.*((-1)^(i-1))/(2.*i-1);
    pi_2(i)=pi_2(i-1)+(factorial(i-1)^2)*(2.^i)/factorial(2*i-1);
end

figure(3)
plot(pi_1)
xlabel('n', 'FontSize', 30)
ylabel('\pi_1', 'FontSize', 30)

figure(4)
plot(pi_2)
xlabel('n', 'FontSize', 30)
ylabel('\pi_2', 'FontSize', 30)

% Problem 10
% a)
% x coordinate where the taylor series is expanded around
x=3.;
% first term
e(1)=1.;
for i=2:15
   e(i)=e(i-1)+(x^(i-1))/factorial(i-1);
   err(i-1)=abs((e(i)-e(i-1))/e(i));
end

% b)
% tolerance value
tol=1E-4;
% error
err(1)=10.^10;
% the first indice
i=1;
 while err(i,1) > tol
   i=i+1;
   e(i)=e(i-1)+(x^(i-1))/factorial(i-1);
   err(i)=abs((e(i)-e(i-1))/e(i));
 end