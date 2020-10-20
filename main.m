clear all
%setenv('MW_MINGW64_LOC','E:\mingw-w64\x86_64-5.3.0-posix-seh-rt_v4-rev0\mingw64');
clc
%mex cec20_func.cpp -DWINDOWS
D=5;%维度
Xmin=-100;
Xmax=100;
pop_size=100;
runs=30;
switch D
    case 5
        EFS=50000;
    case 10
        EFS=1000000;
    case 15
        EFS=3000000;
    case 20
        EFS=10000000;
    otherwise
        a='维度不存在';
        fprintf('%s',a);
end
fhd=str2func('cec20_func');

TbL_SHADE_D5=zeros(16,runs,10);
TbL_SHADE_D5_result=zeros(10,5);%存储每次重复每个函数的最优值
fbest=zeros(10,runs);
pop_anay=zeros(runs,3,10);
TbL_SHADE_anay_result_D5=zeros(10,3);%种群聚类和种群密度结果
for i=1:10%函数序号
    func_num=i;
    for j=1:runs
        i,j,
        [bestvalue,TbL_SHADE_D5(:,j,i),pop_anay(j,:,i)]=TbL_SHADE(fhd,D,pop_size,Xmin,Xmax,EFS,func_num);
        fbest(i,j)=bestvalue;
    end
    TbL_SHADE_D5_result(i,1)=min(fbest(i,:));
    TbL_SHADE_D5_result(i,2)=max(fbest(i,:));
    TbL_SHADE_D5_result(i,3)=median(fbest(i,:));
    TbL_SHADE_D5_result(i,4)=mean(fbest(i,:));
    TbL_SHADE_D5_result(i,5)=std(fbest(i,:));
end
for i=1:10
    temp_1=pop_anay(:,:,i);
    a=sum(temp_1(:,1)==1);%出现聚类的实验次数
    TbL_SHADE_anay_result_D5(i,1)=a;
    if(a>=1)
        temp_2=temp_1(temp_1(:,1)==1,:);
        TbL_SHADE_anay_result_D5(i,2)=mean(temp_2(:,2));
        TbL_SHADE_anay_result_D5(i,3)=mean(temp_2(:,3));
    end
end



