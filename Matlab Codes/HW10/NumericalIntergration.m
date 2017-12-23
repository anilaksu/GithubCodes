% this script takes numerical integral

% Problem 1
% the start point of the integration
a=2.;
% the end point of the integration
b=5.;
% the number intervals
n=1000;
% the step size
dx=(b-a)/n;
% the interval
x=linspace(a,b,n+1);

% a) Trapezoid Method

for i=1:n
    % the value at the left point of the interval
   [ f_l ] = getF( x(i) );
    % the value at the left point of the interval
   [ f_r ] = getF( x(i+1) );
   % the area under each trapezoid
   I_i(i)=0.5*dx*(f_l+f_r);
end

% the total area by Trapezoidal rule
Area_trap=sum(I_i);

% b) Simpson's rule by 3/8

for i=1:n
   % the value at the left point of the interval
   [ f_l ] = getF( x(i) );
   % the first  3/8 
   [ f_f ] = getF( (2*x(i)+x(i+1))/3. );
   % the second 3/8
   [ f_s ] = getF( (x(i)+2.*x(i+1))/3. );
    % the value at the left point of the interval
   [ f_r ] = getF( x(i+1) ); 
   % the area under each trapezoid
   I_i(i)=dx*(f_l+3.*f_f+3.*f_s+f_r)/8;
end

% the total area by Simpson's rule
Area_Simpson=sum(I_i);

% c) Gauss Quadrature

for i=1:n
   % the value at the left point of the interval
   [ f_l ] = getF( x(i) );
   % the first quadrature point
   [ f_f ] = getF( 0.5*(x(i)+x(i+1))-((3/5)^(0.5))*(x(i+1)-x(i))/2. );
   % the second quadrature point
   [ f_s ] = getF( 0.5*(x(i)+x(i+1))+((3/5)^(0.5))*(x(i+1)-x(i))/2.);
   % the area under each trapezoid
   I_i(i)=0.5*dx*(8.*f_l+5.*f_f+5.*f_s)/9;
end

% the total area by Gauss Quadrature
Area_Gauss=sum(I_i);

% Problem 2

% the friction coefficient
r=2.;
% the gravity constan
g=9.81;
% the number of time steps
N=100;
% the time step size
dt=0.1;
% the velocity array
v=zeros(N,1);
% the time array
t=zeros(N,1);
% the velocity at the first time step
v(1)=10.;
% let's start the time integration,
for i=1:(N-1)
   v(i+1)=v(i)+dt*(g-r*v(i));
   t(i+1)=t(i)+dt;
end

figure(1)
plot(t,v)
xlabel('time(sec)')
ylabel('Velocity')

plot(t,v)


