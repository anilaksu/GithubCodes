% this script is developed to solve HW6

%Problem 1

% tolerance
tol=1E-6;
% the radius 
R=3.;
% the initial guess
h_0=1.;
% the volume
V_0=30.;
% the error 
err_1=1.;

while err_1 > tol
% the derivative of the volume function
[ df ] = getDerivativeVolume( h_0,R  ) ;
% the volume function
[ f ] = getVolume( V_0,h_0,R );
% the next guess
h_n=h_0-f/df;
% the absolute relative error
err_1=abs((h_n-h_0)/h_n);
% the next iteration value
h_0=h_n;

end

% Problem 2

% x_0
x_0=0.5;
% the step size
dx=0.33;
% the order of the Taylor series
n=15;
% the sum of the Taylor series
sum_T=0.;

for i=0:n
    % the i-th derivative
    [ f_i ] = getNthDerivative( i,x_0 )
    sum_T=sum_T+f_i*(dx^i)/factorial(i);
end

% Problem 3

A=[0 10 21 56 ; 3 40 55 1; 62 26 18 89; 2 24 5 96];

% two variables
x_1=A(3,4);
x_2=A(4,4);
% c array
c=A(1,:);
% d array
d=A(:,2);
% e array
e=A(3,1:3);
% f array
f=A(1:2,1);
% G matrix
G=A(1:2,1:2);

% Problem
S=[1 1 -1; 1 2 0; -1 0 5];
% the transpose of S matrix
S_t=transpose(S);
T=S-S_t;

Q=[3/7 2/7 6/7; -6/7 3/7 2/7; 2/7 6/7 -3/7];
Q_t=transpose(Q);
Q_i=inv(Q);
L=Q_t-Q_i;

I=eye(3);

% Problem 5
A=[2 -6 1; -3 -1 7; -8 1 -2];
b=[-38; -34; -20];
% the solution
x=inv(A)*b;

% Problem 6
A=[0 10 20 30; 40 50 60 70; 80 90 0 10; 20 30 40 50];

% the sum of all elements of the third row
sum_3=0.;
% the sum of the diagonal element
sum_d=0.;
for j=1:length(A(:,1))
    sum_3=sum_3+A(3,j);
    sum_d=sum_d+A(j,j);
end

% the sum of the all elements
sum_all=0.;
for i=1:length(A(:,1))
    for j=1:length(A(:,1))
    sum_all=sum_all+A(i,j);
    end
end

% Problem 7 
P=[0 0; 2 1; 4 3];
% L1
L_1=((P(1,1)-P(2,1))^2.+(P(1,2)-P(2,2))^2.)^0.5;
% L2
L_2=((P(2,1)-P(3,1))^2.+(P(2,2)-P(3,2))^2.)^0.5;
% L3
L_3=((P(1,1)-P(3,1))^2.+(P(1,2)-P(3,2))^2.)^0.5;
% L_a
L_a=(L_1+L_2+L_3)/2.;
% the angles
a_1=atan(((P(1,2)-P(2,2))/(P(1,1)-P(2,1))));
a_2=atan(((P(3,2)-P(2,2))/(P(3,1)-P(2,1))));
a_3=atan(((P(3,2)-P(1,2))/(P(3,1)-P(1,1))));
% the area 
Area=(L_a*(L_a-L_1)*(L_a-L_2)*(L_a-L_3))^.5;