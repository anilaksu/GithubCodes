% bisection root finding method

format long

% let's generate our grid 
x=linspace(0.4,0.48,100);

% f(x) function evaluated at the grid points
y=tan(pi*x)-x-6;

% let's generate plot 
figure(1)

plot(x,y);
%title('f(x) ') % title
ylabel('f(x)')
% label for y axis
xlabel('x')
% label for x axis

%Number of iterations
N=50;

% let's define xl vector 
xl=zeros(N);

% let's define xl vector
xr=zeros(N);

% first left point
xl(1)=0.4;

%first right point
xr(1)=0.48;

for i=1:N-1
    
   % by shorthand notation functions value at left point
   fl=tan(pi*xl(i))-xl(i)-6;
   
   % by shorthand notation functions value at left point
   fr=tan(pi*xr(i))-xr(i)-6;
   
   % middle point
   xm=0.5*(xl(i)+xr(i));
   
   % function evaluated at that point
   fm=tan(pi*xm)-xm-6
   
    if fl*fm<0 
        xr(i+1)=xm;
    
        xl(i+1)=xl(i);
        
    else
        xl(i+1)=xm;
    
        xr(i+1)=xr(i); 
        
    end    
        
     % this is our cut-off criteria
        if abs(fl-fm)<10^-5.
            
            break
        end
end
    

xl(i)

xr(i)

fl=tan(pi*xr(i))-xr(i)-6

i

    
