% this program is written to process the spectrometer measurements

format long

%% Argon Spectrometer Data 
ArgonData = xlsread('20161110OESdata1.xlsx')

%% kb constant
kb=1.380662*10^(-23);

% wavelength of Argon Gas
lambda=zeros(length(ArgonData(:,1)),1);
% statistical weight
g=zeros(length(ArgonData(:,1)),1);
% Intensity
I=zeros(length(ArgonData(:,1)),1);
% transition probability
A=zeros(length(ArgonData(:,1)),1);
% excitation energy
Ek=zeros(length(ArgonData(:,1)),1);

% let's distribute the read data
lambda=ArgonData(:,1);
g=ArgonData(:,2);
I=ArgonData(:,3);
A=ArgonData(:,4);
Ek=ArgonData(:,5);

% let's define new variables y and x
for i=1:length(ArgonData(:,1))
    % y is defined by the formula below
    y(i)=log((I(i)*lambda(i))/(g(i)*A(i)));
    % x is defined by the formula below
    x(i)=Ek(i)/kb;
end

%% we are trying to fit y=c1+c2x by means of least square
A=ones(length(ArgonData(:,1)),2);
%% we set first column 1 and second colum is equal to x
A(:,2)=x;
% %% let's transform min:sec to sec
% time=zeros(length(ExpData(:,1)),1);
% 
% %% It is incremental time
% for i=1:length(ExpData(:,1))
%     time(i)=ExpData(i,25)-ExpData(1,25)+60*(ExpData(i,24)-ExpData(1,24));
% end
% 
figure(1)
 plot(x,y,'r*')
% hold on
% plot(time,ExpData(:,5),'m')
% 
% 
% title('First-Phase Ua and Ia')
% legend('Ua','Ia')
% % label for y axis
% xlabel('time in sec')