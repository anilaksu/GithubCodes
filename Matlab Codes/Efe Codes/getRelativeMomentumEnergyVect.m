function [ p,T ] = getRelativeMomentumEnergyVect( m,v )
%This function is written to relatistic momentum and energy
% speed of light Note that it is m/s
c= 299*10^6
% since we receive v as v/c let's correct it
v=v*c;
% relativistic momentum
p=zeros(length(v),1)
% relativistic energy
T=zeros(length(v),1)

for i=1:length(v)
    %relativistic momentum
    p(i,1)=(m*v)/sqrt(1-(v(i)/c).^2);
    % relativistic energy
    T(i,1)=(m*c^2)/sqrt(1-(v(i)/c).^2)-m*c^2;
end

end
