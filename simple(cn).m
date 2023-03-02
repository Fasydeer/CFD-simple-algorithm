%写在前面：matlab一般无法吃满cpu，因此运算时间较长，可以通过暂停程序查看参数c来看已经迭代次数

m=101;

n=101;%修改除个位数外的位数来修改网格数的数量

P=zeros(m,n);

PX=P;

PP=P;

U=zeros(m,n+1);

V=zeros(m+1,n);

U(1,:)=1;%U的顶盖

VS=V;

US=U;%US是U星

RE=400;%高雷诺数建议增加网格数量进行运算

x=1/(m-1);

y=1/(n-1);

t=0.005;%△xyt的大小，△t<=△x收敛的可能性大，若e=NAN需将△t调小，例如Re=5000时可以调成0.001

e=1;

e1=1;%e1为d的最大值

e2=1;%e2用于选取各个d的值

e3=0;%e3用来保护e1

c=0;%外部迭代次数

l=0;%内部迭代次数

S1=P;%将U从交错网格上回归到P所在的原网格上

S2=P; %将V从交错网格上回归到P所在的原网格上

q=0;%总的内部循环次数

while  c<30000 %%%重要部分！！重要部分！！重要部分！！重要部分！！重要部分！！重要部分！！c<30000是流场收敛的大体迭代次数，已经可以满足大部分情况，也可以自己修改迭代次数缩短运行时间。若对于精度有极高的要求，可以将c<30000改为e1>0.0001，注意，这将消耗大量的运算时间。

     for a=2:n-1

    for b=2:m

  US(a,b)=U(a,b)-t*((U(a,b+1)^2-U(a,b-1)^2)/(2*x)+(U(a+1,b)*(V(a+1,b-1)+V(a+1,b))-U(a-1,b)*(V(a,b-1)+V(a,b)))/(4*y))+(t/RE)*((U(a,b+1)-2*U(a,b)+U(a,b-1))/(x*x)+(U(a+1,b)-2*U(a,b)+U(a-1,b))/(y*y))-(P(a,b)-P(a,b-1))*(t/x);

  VS(b,a)=V(b,a)-t*((V(b+1,a)^2-V(b-1,a)^2)/(2*y)+(V(b,a+1)*(U(b-1,a+1)+U(b,a+1))-V(b,a-1)*(U(b-1,a)+U(b,a)))/(4*x))+(t/RE)*((V(b,a+1)-2*V(b,a)+V(b,a-1))/(x*x)+(V(b+1,a)-2*V(b,a)+V(b-1,a))/(y*y))-(P(b,a)-P(b-1,a))*(t/y); 

    end

     end   %公式1与公式2

   US(:,1)=-US(:,2);

   US(:,n+1)=-US(:,n);

   VS(1,:)=-VS(2,:);

   VS(m+1,:)=-VS(m,:);

    PP=zeros(m,n);%每次内部迭代都要设置一个新的初始迭代值，于是PP于此更新

    PX=PP;

    q=q+l;%总的内部迭代次数

    l=0;

           while 1  %双while语句将无法实现套用的功能，需换成break语句

       l=l+1; %记录单次内部迭代次数     

     for d=2:n-1

    for f=2:m-1   

   PX(d,f)=1/4*(PP(d+1,f)+PP(d,f+1)+PX(d-1,f)+PX(d,f-1)-x*(US(d,f+1)-US(d,f)+VS(d+1,f)-VS(d,f))/t);%等式右边有PX为高斯-塞德尔迭代

    end

     end  %公式6-103

PX(1,:)=PX(2,:);

PX(:,1)=PX(:,2);

PX(m,:)=PX(m-1,:);

PX(:,n)=PX(:,n-1);

   e=max(max(abs(PP-PX)));%检验内部迭代收敛

   E(q+l,1)=e;%将每个误差存放于一个列向量中方便查看收敛趋势

   PP=PX;

   PX=PP;%高斯塞德尔迭代需要这一步来使第一步迭代准确

       if e<=0.0001

           break;%收敛的精度

       end

          if l>2000

           break;%收敛的迭代次数上限

          end

           end %求公式3的迭代

   P=P+0.1*PX;   %得到当前时刻的P值

P(1,:)=P(2,:);

P(:,1)=P(:,2);

P(m,:)=P(m-1,:);

P(:,n)=P(:,n-1); 

e3=0;

   for g=2:n-1

    for h=2:m-1   

   e2=abs(-US(g,h+1)/x+US(g,h)/x-VS(g+1,h)/x+VS(g,h)/x);

         if e2>e3

             e3=e2;

             e1=e2; %求d的最大值来与精度进行比较

         end

    end

   end  %公式3中d趋近于0

 for a=2:n-1

    for b=2:m   

  U(a,b)=US(a,b)-t/x*(PX(a,b)-PX(a,b-1));

  V(b,a)=VS(b,a)-t/y*(PX(b,a)-PX(b-1,a));

    end

 end %UV的修正方程即公式4

  U(:,1)=-U(:,2);

   U(:,n+1)=-U(:,n);

   V(1,:)=-V(2,:);

   V(m+1,:)=-V(m,:);

     for a=2:n-1

    for b=2:m

  US(a,b)=U(a,b)-t*((U(a,b+1)^2-U(a,b-1)^2)/(2*x)+(U(a+1,b)*(V(a+1,b-1)+V(a+1,b))-U(a-1,b)*(V(a,b-1)+V(a,b)))/(4*y))+(t/RE)*((U(a,b+1)-2*U(a,b)+U(a,b-1))/(x*x)+(U(a+1,b)-2*U(a,b)+U(a-1,b))/(y*y))-(P(a,b)-P(a,b-1))*(t/x);

  VS(b,a)=V(b,a)-t*((V(b+1,a)^2-V(b-1,a)^2)/(2*y)+(V(b,a+1)*(U(b-1,a+1)+U(b,a+1))-V(b,a-1)*(U(b-1,a)+U(b,a)))/(4*x))+(t/RE)*((V(b,a+1)-2*V(b,a)+V(b,a-1))/(x*x)+(V(b+1,a)-2*V(b,a)+V(b-1,a))/(y*y))-(P(b,a)-P(b-1,a))*(t/y); 

    end

     end  %将修正好的UVP再代入修正方程中求得下一时刻的UV

   US(:,1)=-US(:,2);

   US(:,n+1)=-US(:,n);

   VS(1,:)=-VS(2,:);

   VS(m+1,:)=-VS(m,:);

U=US;

V=VS;

c=c+1;%外部迭代次数

end

for a=1:n

    for b=1:m 

        S1(a,b)=(US(a,b)+US(a,b+1))/2;

        S2(b,a)=(VS(b,a)+VS(b+1,a))/2;%将UV从交错网格上重新放到P的网格中

    end

end %将UV化成101*101网格

xx=0:m-1;

yy=0:n-1;

[X,Y]=meshgrid(xx,yy);%生成坐标

figure
pcolor(X,Y,P);%生成压强图

shading interp%输出图像

figure

streamslice(S1,S2);%生成流线图
