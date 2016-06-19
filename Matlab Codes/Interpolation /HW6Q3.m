% PS06 Q3

x=linspace(-1,1,101);

% let's generate the actual function
for i=1:101
    
    f(i)=1./(1.+25.*x(i)^2.);
    
end

n=8;
xdata=linspace(-1,1,n+1);

% let's generate the actual function
for i=1:n+1
    
    fdata(i)=1./(1.+25.*xdata(i)^2.);
    
end


% let's interpolate it 
[fout] = interphomecooked(xdata,fdata,n,x);


figure(4)
plot(x,fout(:,1),'LineWidth',2)
%title('f(x) vs T_N(x) ')
%legend('f(x)','L_5{}_{}_,_3(x)','L_5{}_{}_,_4(x)','L_5{}_{}_,_5(x)',4)
ylabel('P_4{}_0(x)')
% label for y axis
xlabel('x')
% label for x axis


% let's generate the chebyshev points

n=8;

for i=0:n
    
  xcheb(i+1)=cos((2.*i+1.)*pi/(2.*n+2.));
    
  fcheb(i+1)=1./(1.+25.*xcheb(i+1)^2.);
  
end

[foutc] = interphomecooked(xcheb,fcheb,n,x);

figure(5)

plot(x,foutc(:,1),'LineWidth',2)
%title('f(x) vs T_N(x) ')
%legend('f(x)','L_5{}_{}_,_3(x)','L_5{}_{}_,_4(x)','L_5{}_{}_,_5(x)',4)
ylabel('P_4{}_0(x)')
% label for y axis
xlabel('x')
% label for x axis

% let's calculate the error for both chebysev and equidistant points

for i=1:101

    % absolute error for equidistant mesh points
erreq(i)=abs(fout(i,1)-f(i));

    % absolute error for chebysev points
errcheb(i)=abs(foutc(i,1)-f(i));    


end 

figure(6)
plot(x,erreq,'LineWidth',2)
title('Absolute Error for Equally Spaced Interpolation Points ')
%legend('f(x)','L_5{}_{}_,_3(x)','L_5{}_{}_,_4(x)','L_5{}_{}_,_5(x)',4)
ylabel('Absolute Error')
% label for y axis
xlabel('x')

figure(7)
plot(x,errcheb,'LineWidth',2)
title('Absolute Error for Chebyshev Interpolation Points ')
%legend('f(x)','L_5{}_{}_,_3(x)','L_5{}_{}_,_4(x)','L_5{}_{}_,_5(x)',4)
ylabel('Absolute Error')
% label for y axis
xlabel('x')
