%Matlab generally cannot fill up the cpu, so the calculation time is longer, you can check the parameter "c" by pausing the program to see the number of iterations

m=101;

n=101;%Modify the number of hundreds to modify the number of grid numbers

P=zeros(m,n);

PX=P;

PP=P;

U=zeros(m,n+1);

V=zeros(m+1,n);

U(1,:)=1;%lid speed

VS=V;

US=U;%US is U star

RE=400;%For high Reynolds numbers, it is recommended to increase the number of grids for calculation

x=1/(m-1);

y=1/(n-1);

t=0.005;%when △t<=△x,the possibility of convergence is higher，if e=NAN △t need to be decreased，like when Re=5000,△t needs to be reduced to 0.001

e=1;

e1=1;%e1 is the maximum value of d

e2=1;%e2 is used to replace d

e3=0;%e3is used to "protect" e1

c=0;%number of external iterations

l=0;%number of internal iterations

S1=P;%Return U from the staggered grid to the original grid where P is located

S2=P; %Return V from the staggered grid to the original grid where P is located

q=0;%total number of internal iterations

while  c<30000 %%%important part! ! important part! ! important part! ! important part! ! important part! ! important part! ! c<30000 is the general number of iterations for flow field convergence, which can satisfy most situations, and you can also modify the number of iterations to shorten the running time. 
%If you have extremely high requirements for precision, you can change c<30000 to e1>0.0001. Note that this will consume a lot of computing time.

     for a=2:n-1

    for b=2:m

  US(a,b)=U(a,b)-t*((U(a,b+1)^2-U(a,b-1)^2)/(2*x)+(U(a+1,b)*(V(a+1,b-1)+V(a+1,b))-U(a-1,b)*(V(a,b-1)+V(a,b)))/(4*y))+(t/RE)*((U(a,b+1)-2*U(a,b)+U(a,b-1))/(x*x)+(U(a+1,b)-2*U(a,b)+U(a-1,b))/(y*y))-(P(a,b)-P(a,b-1))*(t/x);

  VS(b,a)=V(b,a)-t*((V(b+1,a)^2-V(b-1,a)^2)/(2*y)+(V(b,a+1)*(U(b-1,a+1)+U(b,a+1))-V(b,a-1)*(U(b-1,a)+U(b,a)))/(4*x))+(t/RE)*((V(b,a+1)-2*V(b,a)+V(b,a-1))/(x*x)+(V(b+1,a)-2*V(b,a)+V(b-1,a))/(y*y))-(P(b,a)-P(b-1,a))*(t/y); 

    end

     end   %Formula 1 and Formula 2

   US(:,1)=-US(:,2);

   US(:,n+1)=-US(:,n);

   VS(1,:)=-VS(2,:);

   VS(m+1,:)=-VS(m,:);

    PP=zeros(m,n);%A new initial iteration value is set for each internal iteration, so PP is updated here

    PX=PP;

    q=q+l;

    l=0;

           while 1  %Double while statement will not be able to implement the applied function, it needs to be replaced with a break statement

       l=l+1; %Record the number of  internal iterations in single external iteration
     for d=2:n-1

    for f=2:m-1   

   PX(d,f)=1/4*(PP(d+1,f)+PP(d,f+1)+PX(d-1,f)+PX(d,f-1)-x*(US(d,f+1)-US(d,f)+VS(d+1,f)-VS(d,f))/t);%PX on the right side of the equation is Gauss-Seidel iteration

    end

     end  

PX(1,:)=PX(2,:);

PX(:,1)=PX(:,2);

PX(m,:)=PX(m-1,:);

PX(:,n)=PX(:,n-1);

   e=max(max(abs(PP-PX)));%Check inner iteration convergence

   E(q+l,1)=e;%Store each error in a column vector for easy viewing of convergence trends

   PP=PX;

   PX=PP;%Gauss-Seidel iteration needs this step to make the first iteration accurate

       if e<=0.0001

           break;%Convergence precision
       end

          if l>2000

           break;%Maximum number of iterations to converge

          end

           end %Find iterations of Equation 3

   P=P+0.1*PX;   %Get the P value at the current moment
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

             e1=e2; %Find the maximum value of d to compare with precision

         end

    end

   end  %d tends to 0 in formula 3

 for a=2:n-1

    for b=2:m   

  U(a,b)=US(a,b)-t/x*(PX(a,b)-PX(a,b-1));

  V(b,a)=VS(b,a)-t/y*(PX(b,a)-PX(b-1,a));

    end

 end %The UV correction equation is formula 4

  U(:,1)=-U(:,2);

   U(:,n+1)=-U(:,n);

   V(1,:)=-V(2,:);

   V(m+1,:)=-V(m,:);

     for a=2:n-1

    for b=2:m

  US(a,b)=U(a,b)-t*((U(a,b+1)^2-U(a,b-1)^2)/(2*x)+(U(a+1,b)*(V(a+1,b-1)+V(a+1,b))-U(a-1,b)*(V(a,b-1)+V(a,b)))/(4*y))+(t/RE)*((U(a,b+1)-2*U(a,b)+U(a,b-1))/(x*x)+(U(a+1,b)-2*U(a,b)+U(a-1,b))/(y*y))-(P(a,b)-P(a,b-1))*(t/x);

  VS(b,a)=V(b,a)-t*((V(b+1,a)^2-V(b-1,a)^2)/(2*y)+(V(b,a+1)*(U(b-1,a+1)+U(b,a+1))-V(b,a-1)*(U(b-1,a)+U(b,a)))/(4*x))+(t/RE)*((V(b,a+1)-2*V(b,a)+V(b,a-1))/(x*x)+(V(b+1,a)-2*V(b,a)+V(b-1,a))/(y*y))-(P(b,a)-P(b-1,a))*(t/y); 

    end

     end  %Substitute the corrected UVP into the correction equation to obtain the UV at the next moment

   US(:,1)=-US(:,2);

   US(:,n+1)=-US(:,n);

   VS(1,:)=-VS(2,:);

   VS(m+1,:)=-VS(m,:);

U=US;

V=VS;

c=c+1;%number of external iterations

end

for a=1:n

    for b=1:m 

        S1(a,b)=(US(a,b)+US(a,b+1))/2;

        S2(b,a)=(VS(b,a)+VS(b+1,a))/2;%Reposition the UVs from the staggered mesh to P's mesh

    end

end %Convert UV into 101*101 grid

xx=0:m-1;

yy=0:n-1;

[X,Y]=meshgrid(xx,yy);%generate coordinates

figure
pcolor(X,Y,P);%Generate a pressure map

shading interp%output image

figure

streamslice(S1,S2);%Generate streamlines
