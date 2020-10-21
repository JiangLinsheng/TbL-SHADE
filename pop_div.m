% Copyright (c) 2020, Linsheng Jiang
% All rights reserved.
%Calculate population density
function PD=pop_div(pop,Np,D)
    xj=mean(pop,1);%the average value of each dimension
    for j=1:D
        pop(:,j)=pop(:,j)-xj(j);
    end
    pop=power(pop,2);
    PD=sqrt(sum(pop(:))/Np);
end
