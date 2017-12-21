% this script is developed to solve HW6

%Problem 1

% the number of iterations
N=5;
% the x vector
x=zeros(2,N+1);
% the inital guess
x(1,1)=-0.1;
x(2,1)=1.0;
% f vector
f=zeros(2,N+1);
% the Jacobian
J=zeros(2,2);

for i=1:N
    % the evaluation of function at i-th iteration
    f(1,i)=10.-x(1,i)^2.-x(2,i)^2.;
    f(2,i)=1.-2.*x(1,i)-x(2,i);
    % the jacobian matrix
    J(1,1)=-2.*x(1,i);
    J(1,2)=-2.*x(2,i);
    J(2,1)=-2.;
    J(2,2)=-1.;
    % x vector at next vector
    x(:,i+1)=x(:,i)-inv(J)*f(:,i);
end

% Problem 2
P=[47 89; 50 0; 13.5 16; 0 50; 0 127.5; 20 147.5; 79 111; 100 50; 84 13];
% Area vector
Area=zeros(8,1);
for i=1:8
    if(i<8)
        % L1
        L_1=((P(1,1)-P(i+1,1))^2.+(P(1,2)-P(i+1,2))^2.)^0.5;
        % L2
        L_2=((P(i+1,1)-P(i+2,1))^2.+(P(i+1,2)-P(i+2,2))^2.)^0.5;
        % L3
        L_3=((P(1,1)-P(i+2,1))^2.+(P(1,2)-P(i+2,2))^2.)^0.5;
        % L_a
        L_a=(L_1+L_2+L_3)/2.;
    else
         % L1
        L_1=((P(1,1)-P(9,1))^2.+(P(1,2)-P(9,2))^2.)^0.5;
        % L2
        L_2=((P(9,1)-P(2,1))^2.+(P(9,2)-P(2,2))^2.)^0.5;
        % L3
        L_3=((P(2,1)-P(1,1))^2.+(P(2,2)-P(1,2))^2.)^0.5;
        % L_a
        L_a=(L_1+L_2+L_3)/2.;
    end
% the area 
Area(i)=(L_a*(L_a-L_1)*(L_a-L_2)*(L_a-L_3))^.5; 
end

% Problem 3

% x vector
x=linspace(1,6,11);

for i=1:10
   f_i = getF( x(i) );
   f_p = getF( x(i+1) );
   I_i(i)=0.25*(f_i+f_p); 
end

% the area
Area_t=sum(I_i);

% Problem 4
COOR=dlmread('LSR.txt');