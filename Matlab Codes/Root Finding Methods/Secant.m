% Secant root finding method
% created by Anil Aksu
% I don't claim copyright you can use it for any purpose you need

format long

%Number of iterations
N=25;

% iteration point 
x=zeros(N);

% function at that point
f=zeros(N);

%relative error vector
Err=zeros(N-2);

%initial guess
x(1)=0.54;

% the second initial point to form secant line
x(2)=0.48;


for i=2:N-1
    
   % function evaluated at first point
   f(i)=tan(pi*x(i))-x(i)-6;
   
   %first derivative evaluated at that poin
   f(i-1)=tan(pi*x(i-1))-x(i-1)-6;
   
   % next iteration point
   x(i+1)=x(i)-f(i)*(x(i)-x(i-1))/(f(i)-f(i-1));
   
  % relative error
   Err(i)=abs(x(i+1)-x(i));
        
     % this is our cut-off criteria
        if Err(i)<10^-6.
            
            break
        end
end
    
loglog(Err);
%title
title('Relative Error for Secant Method')

ylabel('Relative Error')
% label for y axis
xlabel('Number of iteration')
% label for x axis

i

fl=tan(pi*x(i))-x(i)-6



    
