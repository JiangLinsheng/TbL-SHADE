%������Ⱥ�ܶ�
function PD=pop_div(pop,Np,D)
    xj=mean(pop,1);%ÿһά��ƽ��ֵ
    for j=1:D
        pop(:,j)=pop(:,j)-xj(j);
    end
    pop=power(pop,2);
    PD=sqrt(sum(pop(:))/Np);
end