function [ df ] = getDerivativeVolume( h,R  )
%this function computes the derivatice of the volume
df=pi*h^2.-2.*pi*h*R;

end

