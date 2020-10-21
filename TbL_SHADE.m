% Copyright (c) 2020, Linsheng Jiang
% All rights reserved.
%L-SHADE using turning-based mutation
function [g_best,error_value_1,pop_anay]=TbL_SHADE(fhd,D,pop_size,Xmin,Xmax,EFS,varargin)
format long;
rand('state',sum(100*clock));
Np = pop_size;%Population size
Np_init=pop_size;
Np_min=4;
xmin = Xmin;%lower bound
xmax = Xmax;%upper bound
value = zeros(1,Np);%fitness value
value_1=zeros(1,Np);%Store the fitness value of the previous generation, because sort() changes the value array
pop_next_1=zeros(Np,D);%Population after mutation and crossover
pop_next=zeros(Np,D);%Next generation population after selection
A=zeros(Np,D);
A_pos=0;
g=1;%the number of generation
efs=0;%the number of the calculation of the fitness function
min_error=1000;
error_value_1=zeros(16,1,1);%Statistical experiment results
memory_pos=1;
k=0;
optimum=[100,1100,700,1900,1700,1600,2100,2200,2400,2500];%The optimal value of 10 test functions
%the parameters of population diversity analysis
pop_anay=zeros(1,3,1);
epsilon=2;
MinPts=4;
function y = f(pos)
    y=feval(fhd,pos',varargin{:});
    error_y=y-optimum(varargin{:});
    if error_y<power(10,-8)
        error_y=0;
        min_error=0;
    end
    efs=efs+1;%the number of the calculation of the fitness function self adding
    if efs<=EFS
        spefs=round(power(D,(k/5-3))*EFS);
        if efs == spefs
            k=k+1;
            error_value_1(k,1,1)=error_y;
        end
    end
end

%the parameters of turning-based mutation
distance_init=0;%the initial outer radius
distance_min=0;%the inner radius
for di=1:D
    distance_init=distance_init+power(100,2);
    distance_min=distance_min+power(5,2);
end
distance_init=sqrt(distance_init);
distance_min=sqrt(distance_min);

pop=rand(Np,D)*(xmax-xmin)+xmin;% initialize the population
for m=1:Np%Calculate the fitness value
    value(m)=f(pop(m,:));
end
g_best=min(value);%Global optimal value
picture_data(g)=g_best;%Data for drawing
value_1=value;
Mf=zeros(1,Np);%the historical memory of scaling factor F
Mcr=zeros(1,Np);%the historical memory of crossover rate CR
Mf=Mf+0.5;
Mcr=Mcr+0.5;
while (efs<EFS && min_error>=power(10,-8))
    if(pop_anay(1,1,1)~=1)
        [~, isnoise]=DBSCAN(pop,epsilon,MinPts);
        if(any(isnoise==false))
           pop_anay(1,1,1)=1;
           pop_anay(1,2,1)=g;
           pop_anay(1,3,1)=pop_div(pop,Np,D);
        end
    end
    
    Sf=zeros(1,Np);%Hold the better scaling factor F in each generation
    Scr=zeros(1,Np);%Hold the better crossover rate CR in each generation
    dif_fitness=zeros(1,Np);%the difference of fitness
    dif_fitness_sum=0;%the sum of the difference of fitness
    dif_sum=0;%the number of difference
    [~,b]=sort(value,2);%sort
    distance=round(((distance_min-distance_init)/EFS)*efs+distance_init);
    for i=1:Np
        
        %%%%%%%%%%%%%%%%%%%%%%%%----mutation----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cri=randi(Np);
        fi=cri;
        CRi=normrnd(Mcr(1,cri), 0.1);%Normal distribution
        if(CRi>=1)
            CRi=1;
        elseif(CRi<=0)
            CRi=0;
        end
        Fi= Mf(1,fi) + 0.1 * tan(pi * (rand - 0.5));%Cauchy distribution
        while(Fi<=0)
            Fi= Mf(1,fi) + 0.1 * tan(pi * (rand - 0.5));
        end
        if(Fi>=1)
            Fi=1;
        end
        
        ps=rand*(0.2-2/Np)+2/Np;
        Pbest=b(randi(round(Np*ps)));%the individual selected from the best NP × ps (ps ∈ [0,1]) individuals randomly in the current population 
         dx = randperm(Np);          
         r1=dx(1);%The other two individuals
         if r1 == i
            r1= dx(2);
         end
         pop_1=[pop;A];
         if(A_pos<=Np)
            r2g=round(randi(Np+A_pos));
         else
             r2g=round(randi(Np*2));
         end
         while(r2g==i||r2g==r1)
             if(A_pos<=Np)
                r2g=round(randi(Np+A_pos));
             else
                r2g=round(randi(Np*2));
             end
         end
         
         par_distance=0;
         for dj=1:D
             par_distance=par_distance+power((pop(i,dj)-pop(Pbest,dj)),2);
         end
         par_distance=sqrt(par_distance);
         temp_de=Fi*(pop(Pbest,:)-pop(i,:))+Fi*(pop(r1,:)-pop_1(r2g,:));%Difference vector
         if par_distance<distance && par_distance>distance_min
             temp_de=-temp_de;
             temp_d=randi(D);%the number of randomly selected dimensions
             rand_d=randperm(D);%the randomly disordered dimension index array
             for d=1:temp_d
                 temp_de(rand_d(d))=rand*(xmax-xmin)+xmin;
             end
         end
         son=pop(i,:)+temp_de;%turning based mutation
         
          %%%%%%%%%%%%%%%%%%----boundary-handling and crossover---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         jrand=randi(D);
         for j = 1: D
            while(son(1,j)<xmin || son(1,j)>xmax)
                 if(son(1,j)<xmin)
                     son(1,j) =(pop(i,j)+xmin)/2;
                 end
                 if(son(1,j)>xmax)
                     son(1,j) =(pop(i,j)+xmax)/2;
                 end
            end
                 if jrand==j||rand<=CRi
                    pop_next_1(i,j)=son(1,j);
                 else
                    pop_next_1(i,j)=pop(i,j);
                 end
         end
          %%%%%%%%%%%%%%%%%%----selection---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         old_fit=value_1(i);
         new_fit=f(pop_next_1(i,:));
         if old_fit<=new_fit
            pop_next(i,:)=pop(i,:);
            value(i)=old_fit;
         else
            dif_sum=dif_sum+1;
            pop_next(i,:)=pop_next_1(i,:);
            value(i)=new_fit;
            A_pos=A_pos+1;
            if A_pos<=Np
                A(A_pos,:)=pop(i,:);
            else
                A(randi(Np),:)=pop(i,:);
            end
            Scr(1,dif_sum)=CRi;%the better crossover rate CR
            Sf(1,dif_sum)=Fi;%the better scaling factor F
            dif_fitness(1,dif_sum)=old_fit-new_fit;%the difference of fitness
         end
         value_1(i) = value(i);
    end
    temp_g_best=min(value_1);
    if temp_g_best<g_best
        g_best=temp_g_best;
    end
    %update the historical memories of scaling factor F and crossover rate CR
    if(dif_sum>0)
        dif_fitness_sum=sum(dif_fitness);
        Mcr(1,memory_pos)=0;
        Mf(1,memory_pos)=0;
        temp_sum_sf=0;
        for i=1:dif_sum
            weight=dif_fitness(1,i)/dif_fitness_sum;
            Mcr(1,memory_pos)=Mcr(1,memory_pos)+weight*Scr(1,i);
            Mf(1,memory_pos)=Mf(1,memory_pos)+weight*Sf(1,i)*Sf(1,i);
            temp_sum_sf=temp_sum_sf+weight*Sf(1,i);
        end
        Mf(1,memory_pos)=Mf(1,memory_pos)/temp_sum_sf;
        memory_pos=memory_pos+1;
        if memory_pos>Np_init
            memory_pos=1;
        end
    end
    pop=pop_next;
    %Linear Decrease in Population Size
    Np_next=round((Np_min-Np_init)/EFS*efs+Np_init);
    if(Np_next<Np)
        while(Np>Np_next)
            [~,max_index]=max(value);%sort
            pop(max_index,:)=[];
            A(randi(Np),:)=[];
            value(:,max_index)=[];
            value_1(:,max_index)=[];
            pop_next_1(max_index,:)=[];
            pop_next(max_index,:)=[];
            Np=Np-1;
        end
    end
    g=g+1;
    picture_data(g)=g_best;
end

end
