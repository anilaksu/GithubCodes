function [ p,T ] = getRelativeMomentumEnergy( m,v )
%This function is written to relatistic momentum and energy
% speed of light Note that it is m/s
c= 299*10^6
% since we receive v as v/c let's correct it
v=v*c;
%relativistic momentum
p=(m*v)/sqrt(1-(v/c).^2);
% relativistic energy
T=(m*c^2)/sqrt(1-(v/c).^2)-m*c^2;
end

