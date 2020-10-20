%计算种群密度
function PD=pop_div(pop,Np,D)
    xj=mean(pop,1);%每一维的平均值
    for j=1:D
        pop(:,j)=pop(:,j)-xj(j);
    end
    pop=power(pop,2);
    PD=sqrt(sum(pop(:))/Np);
end