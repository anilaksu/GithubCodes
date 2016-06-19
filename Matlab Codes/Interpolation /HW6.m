% let's generate our mesh points

x=linspace(1,5,100);

% Taylor interpolation function

T=zeros(6,100);

for i=1:100
   % let's calculate function at grid points
   f(i)=1/x(i);
    
   % now let's calculate taylor interpolaiton function
   
   % let's calculate taylor expansion term by term 
   for j=1:6
       
      T(j,i)=((x(i)-1.)^(j-1))*(cos(pi*(j-1)))*(1)^(-1.*j)/factorial(j-1);
       
   end
   
   % now let's add them to generate our interpolation functions
   for j=1:5
       
       T(j+1,i)=T(j+1,i)+T(j,i);
       
   end
   
end



figure(1)

plot(x,f,'LineWidth',2);
hold on
plot(x,T(1,:),'r','LineWidth',2);
hold on
plot(x,T(2,:),'g','LineWidth',2);
hold on
plot(x,T(3,:),'m','LineWidth',2);
grid on
%title
title('f(x) vs T_N(x) ')
legend('f(x)','T_0(x)','T_1(x)','T_2(x)',4)
ylabel('f(x)')
% label for y axis
xlabel('x')
% label for x axis

figure(2)

plot(x,f,'LineWidth',2);
hold on
plot(x,T(4,:),'r','LineWidth',2);
hold on
plot(x,T(5,:),'g','LineWidth',2);
hold on
plot(x,T(6,:),'m','LineWidth',2);
grid on
%title
title('f(x) vs T_N(x) ')
legend('f(x)','T_3(x)','T_4(x)','T_5(x)',4)
ylabel('f(x)')
% label for y axis
xlabel('x')
% label for x axis

% let's calculate percent true error at 3
err=zeros(1,6);
for j=1:6
   
    err(1,j)=abs((f(50)-T(j,50))/f(50))*100.;
    
end

err


% now let's generate interpolation by lagrange interpolants

% interpolant points
xi=[2 2.5 4];

for i=1:100

    
    % now let's calculate the interpolant value
    
    % first lagrange interpolant
    L(1,i)=(x(i)-xi(1))*(x(i)-xi(1))/((xi(1)-xi(2))*(xi(1)-xi(3)));
    
     % second lagrange interpolant
    L(2,i)=(x(i)-xi(2))*(x(i)-xi(2))/((xi(2)-xi(1))*(xi(2)-xi(3)));
    
    % third lagrange interpolant
    L(3,i)=(x(i)-xi(3))*(x(i)-xi(3))/((xi(3)-xi(1))*(xi(3)-xi(2)));
    
    
    % now let's generate the full interpolation polynomial
    P2(i)=L(1,i)*(1/xi(1))+L(2,i)*(1/xi(2))+L(3,i)*(1/xi(3));
    
    
    
end

figure(3)

plot(x,f,'LineWidth',2);
hold on
plot(x,P2,'r','LineWidth',2);
hold on
plot(x,T(3,:),'m','LineWidth',2);
grid on
%title
title('f(x) vs P_2(x) vs T_2(x)')
legend('f(x)','P_2(x)','T_2(x)',3)
ylabel('f(x)')
% label for y axis
xlabel('x')
% label for x axis

errL=abs((f(50)-P2(50))/f(50))*100.
