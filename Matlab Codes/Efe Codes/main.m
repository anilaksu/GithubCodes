% this main.m is generated to call getRelativeMomentumEnergy function
format long
clear all

%Note that we generate ratio v/c here
v=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99,0.999]

% let's pick up arbitrary mass
m=0.1

% relativistic momentum
p=zeros(length(v),1)
T=zeros(length(v),1)
% let's get relativistic momentum and energy
for i=1:length(v)
 [ p(i,1),T(i,1) ] = getRelativeMomentumEnergy( m,v(i) )
end
figure(1)  
plot(v,p(:,1),'r');
title('Relativistic Momentum')

ylabel('Relativistic Momentum')
% label for y axis
xlabel('Velocity v/c')
% label for x axis

figure(2)
plot(v,T(:,1));
title('Relativistic Energy')

ylabel('Relativistic Energy')
% label for y axis
xlabel('Velocity v/c')
% label for x axis