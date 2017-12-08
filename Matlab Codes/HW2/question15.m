function [ d ] = question15( c )
% this function generated d matrix for question 15

for i=1:length(c)
    for j=1:length(c)
        d(i,j)=c(j)+(i-1);
    end
end

end

