% this program is written to process the experimental data of the
% plasmatron

% this data includes time series of currents and voltage differene between
% three phases 

format long

%% sample data 
ExpData=textread('202163.txt')

%% ExpData(:,1)=time 
%% ExpData(:,2)=Ua 
%% ExpData(:,3)=Ub 
%% ExpData(:,4)=Uc 
%% ExpData(:,5)=Ia 
%% ExpData(:,6)=Ib 
%% ExpData(:,7)=Ic
%% ExpData(:,8)=PLin
%% ExpData(:,9)=G water
%% ExpData(:,10)=Pcol
%% ExpData(:,24)=Minutes
%% ExpData(:,25)=Second

%% let's transform min:sec to sec
time=zeros(length(ExpData(:,1)),1);

%% It is incremental time
for i=1:length(ExpData(:,1))
    time(i)=ExpData(i,25)-ExpData(1,25)+60*(ExpData(i,24)-ExpData(1,24));
end

figure(1)
plot(time,ExpData(:,2),'r')
hold on
plot(time,ExpData(:,5),'m')


title('First-Phase Ua and Ia')
legend('Ua','Ia')
% label for y axis
xlabel('time in sec')