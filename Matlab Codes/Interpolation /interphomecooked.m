function [fout] = interphomecooked(xdata,f,n,xout)

% let's calculate the resultant interpolation

% the number of interpolation points
M=length(xout);

% fout resultant interpolation
fout=zeros(M);
% lagrange interpolants
L=ones(n+1,M);

for i=1:M
    
   for j=1:n+1
       
       % let's generate each lagrange interpolant seperately
       for k=1:n+1
          
           if k==j
         L(j,i)=  L(j,i)
           else
              
         L(j,i)=L(j,i)*(xout(i)-xdata(k))/(xdata(j)-xdata(k));      
           end
           
       end
       
   end 
    
   % now let's calculate the resulatant interpolation function 
   
   for j=1:n+1
       
      fout(i)= fout(i)+f(j)*L(j,i);
       
   end
   
end


figure(1)
%plot(xdata,f,'*');
hold on
plot(xout,L(1,:),'r','LineWidth',2);
hold on
plot(xout,L(2,:),'g','LineWidth',2);
hold on
plot(xout,L(3,:),'m','LineWidth',2);
grid on
%title
%title('f(x) vs L_5{}_{}_,_N(x) ')
legend('f(x)','L_5{}_{}_,_0(x)','L_5{}_{}_,_1(x)','L_5{}_{}_,_2(x)',4)
ylabel('f(x)')
% label for y axis
xlabel('x')
% label for x axis

figure(2)

%plot(xdata,f,'*');
hold on
plot(xout,L(4,:),'r','LineWidth',2);
hold on
plot(xout,L(5,:),'g','LineWidth',2);
hold on
plot(xout,L(6,:),'m','LineWidth',2);
grid on
%title
%title('f(x) vs T_N(x) ')
legend('f(x)','L_5{}_{}_,_3(x)','L_5{}_{}_,_4(x)','L_5{}_{}_,_5(x)',4)
ylabel('f(x)')
% label for y axis
xlabel('x')
% label for x axis