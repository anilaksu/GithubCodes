% Newton-Raphson root finding method
% created by Anil Aksu
% I don't claim copyright you can use it for any purpose you need

format long

%Number of iterations
N=25;

% iteration point 
x=zeros(N);

% function at that point
f=zeros(N);

% first derivative at that point
fp=zeros(N);

%relative error vector
Err=zeros(N-1);

%initial guess
x(1)=0.48;

% the root itself
xr=0.4510472588;

for i=1:N-1
    
   % function evaluated at that point
   f(i)=tan(pi*x(i))-x(i)-6;
   
   %first derivative evaluated at that poin
   fp(i)=pi/(cos(pi*x(i))^2)-1;
   
   % next iteration point
   x(i+1)=x(i)-f(i)/fp(i);
   
  % relative error
   Err(i)=abs(x(i+1)-x(i));
        
    % absolute error
   ErrAbs(i)=abs(xr-x(i+1));
   
   % asymptotic error constant
   AsErrCon(i)=abs(xr-x(i+1))/((xr-x(i))^2);
   
     % this is our cut-off criteria
       if Err(i)<10^-9
  
                 break
        end
        
  
        
end
 figure(1)   
loglog(ErrAbs(1:i));
%title
title('Absolute Error for Newton-Raphson Method')

ylabel('Absolute Error')
% label for y axis
xlabel('Number of iteration')
% label for x axis


figure(2)   
loglog(Err(1:i));
%title
title('Relative Error for Newton-Raphson Method')

ylabel('Relative Error')
% label for y axis
xlabel('Number of iteration')
% label for x axis

figure(3)   
plot(AsErrCon(1:i));
%title
title('Asymptotic Error Constant evaluated at each iteration ')

ylabel('Asymptotic Error Constant')
% label for y axis
xlabel('Number of iteration')
% label for x axis

i

fl=tan(pi*x(i))-x(i)-6

AsErrCon(1:i)

% Asymptotic 
AsErr=(pi^2)*sin(pi*xr)/(pi*cos(pi*xr)-(cos(pi*xr)^3))
    
