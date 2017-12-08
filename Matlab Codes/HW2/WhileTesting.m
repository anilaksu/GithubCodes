% this script is developed to test our while skills

vector = [2 4 -5 -3 1 2 -10 2];
% let's define the first indice
i=1;
while i<=length(vector)   
    if(vector(i)<0)
        i=i+1;
        continue
    else
        vector(i)=vector(i)*2;
    end
    i=i+1;
end

% the spatial variables
x=2;
y=3;
z=5;
% let's calculate the radial distance
[ r ] = functionOne( x,y,z );

