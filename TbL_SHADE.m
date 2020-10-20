%L-SHADE based on turning mutations
function [g_best,error_value_1,pop_anay]=TbL_SHADE(fhd,D,pop_size,Xmin,Xmax,EFS,varargin)
format long;
rand('state',sum(100*clock));
Np = pop_size;%��Ⱥ��С
Np_init=pop_size;
Np_min=4;
xmin = Xmin;%��С�½�
xmax = Xmax;%����Ͻ�
value = zeros(1,Np);%��Ӧ��ֵ
value_1=zeros(1,Np);%�洢��һ������Ӧ��ֵ����Ϊsort����ʹvalue�ı���
pop_next_1=zeros(Np,D);%���콻������Ⱥ
pop_next=zeros(Np,D);%ѡ������һ����Ⱥ
A=zeros(Np,D);
A_pos=0;
g=1;%��������
efs=0;%��Ӧ����������
min_error=1000;
error_value_1=zeros(16,1,1);%ͳ��ʵ����
memory_pos=1;
k=0;%ͳ��ʵ������
optimum=[100,1100,700,1900,1700,1600,2100,2200,2400,2500];%ͳ��ʵ������
%��Ⱥ�����Է�����
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
    efs=efs+1;%��Ӧ������������һ
    if efs<=EFS
        spefs=round(power(D,(k/5-3))*EFS);
        if efs == spefs
            k=k+1;
            error_value_1(k,1,1)=error_y;
        end
    end
end

%�������SHADE����
distance_init=0;%��ʼ�������
distance_min=0;%��С�������
for di=1:D
    distance_init=distance_init+power(100,2);
    distance_min=distance_min+power(5,2);
end
distance_init=sqrt(distance_init);
distance_min=sqrt(distance_min);

pop=rand(Np,D)*(xmax-xmin)+xmin;%��ʼ����Ⱥ
for m=1:Np%����Ӧ��ֵ
    value(m)=f(pop(m,:));
end
g_best=min(value);%ȫ������ֵ
picture_data(g)=g_best;
value_1=value;
Mf=zeros(1,Np);%�������ӣ���������)��ʷ����
Mcr=zeros(1,Np);% ����������ʷ����
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
    
    Sf=zeros(1,Np);%�ݴ�ÿ�ε����нϺõ���������
    Scr=zeros(1,Np);%�ݴ�ÿ�ε����нϺõĽ�������
    dif_fitness=zeros(1,Np);%��Ӧ�Ȳ�ֵ
    dif_fitness_sum=0;%��Ӧ�Ȳ�ֵ֮��
    dif_sum=0;%��ֵ����
    [~,b]=sort(value,2);%����
    distance=round(((distance_min-distance_init)/EFS)*efs+distance_init);
    for i=1:Np
        
        %%%%%%%%%%%%%%%%%%%%%%%%----�������----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cri=randi(Np);
        fi=cri;
        CRi=normrnd(Mcr(1,cri), 0.1);%��̬�ֲ�
        if(CRi>=1)
            CRi=1;
        elseif(CRi<=0)
            CRi=0;
        end
        Fi= Mf(1,fi) + 0.1 * tan(pi * (rand - 0.5));%�����ֲ�
        while(Fi<=0)
            Fi= Mf(1,fi) + 0.1 * tan(pi * (rand - 0.5));
        end
        if(Fi>=1)
            Fi=1;
        end
        
        ps=rand*(0.2-2/Np)+2/Np;
        Pbest=b(randi(round(Np*ps)));%ȡ�Ϻý��ǰ%ps�е�һ����
         dx = randperm(Np);          
         r1=dx(1);%������������
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
         temp_de=Fi*(pop(Pbest,:)-pop(i,:))+Fi*(pop(r1,:)-pop_1(r2g,:));%�������
         if par_distance<distance && par_distance>distance_min
             temp_de=-temp_de;
             temp_d=randi(D);%���ѡȡ��ά����
             rand_d=randperm(D);%���ѡȡ��ά��
             for d=1:temp_d
                 temp_de(rand_d(d))=rand*(xmax-xmin)+xmin;%�������
             end
         end
         son=pop(i,:)+temp_de;%����
         
         jrand=randi(D);
         for j = 1: D%�߽紦���������
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
          %%%%%%%%%%%%%%%%%%----ѡ�����---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            Scr(1,dif_sum)=CRi;%�ϺõĽ�������
            Sf(1,dif_sum)=Fi;%�Ϻõı�������
            dif_fitness(1,dif_sum)=old_fit-new_fit;%��Ӧ�Ȳ�ֵ
         end
         value_1(i) = value(i);
    end
    temp_g_best=min(value_1);
    if temp_g_best<g_best
        g_best=temp_g_best;
    end
    %�����������Ӻͽ�������ʷ����
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
    Np_next=round((Np_min-Np_init)/EFS*efs+Np_init);
    if(Np_next<Np)
        while(Np>Np_next)
            [~,max_index]=max(value);%����
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