% this script adds 3 to the elements of array if their value is greater than 3

clear all 
format long 

a =[-5 7 1 0 3 0 5 -10];

% let's check if an element is greater than 3

for i=1:length(a)
    % let's add 3 to those elements
   if(a(i)>3)
       c(i)=a(i)+3;
   else
       c(i)=a(i);
   end
end

% summation of elements of sum_a
sum_a=0;
% let's elements of a vector
for i=1:length(a)    
   sum_a=sum_a+a(i); 
end

% t variable
t=linspace(0,2*pi,201);
% let's try our function
[ f ] = FunctionFor( t );

figure(1)
plot(t,f)
xlabel('t')
ylabel('f(t)')

[ g ] = FunctionWhile( t );

figure(2)
plot(t,g)
xlabel('t')
ylabel('g(t)')


h=t.*exp(t)-t.^2;
 
figure(3)
plot(t,h,'r')
hold on
plot(t,g)
hold on
plot(t,f,'g')

xlabel('t')
%ylabel('h(t)')
legend('h(t)','g(t)','f(t)')

figure(4)
semilogy(t,abs(h),'r')
hold on
semilogy(t,abs(g))
hold on
semilogy(t,abs(f),'g')

xlabel('t')
%ylabel('h(t)')
legend('h(t)','g(t)','f(t)')
