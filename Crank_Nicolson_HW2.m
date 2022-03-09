clc 
clear
close all
%input Datas
k=input('k= ');
h=input('h= ');
C=input('C= ');
A=input('A= ');
P=input('P= ');
T_0=input('T_0= ');
T_inf=input('T_inf= ');
L=input('L= ');
ro=input('density= ');
t=input('time= ');
...Number of pieces created on the length L
J=300;
dx=L/J;
...Time step
ta=0.0005;
t=0:ta:t;
...Calculate alfa and beta
alfa=k/(ro*C);
beta=(h*P)/(A*ro*C);
...Calculation of coefficients a, b, d
a=(alfa*ta/(2*dx^2));
b=((alfa*ta/dx^2))+1;
c=((-alfa*ta/dx^2)+1-(beta*ta));
d=beta.*ta.*T_inf;
...Create a coefficient matrix fot T_n+1
s=zeros(J-1);
for i=1:J-1
    for j=1:J-1
        if i==j
            if i==J-1
               s(i,j)=(b/a)-(1/(1+(h.*dx/k))); 
            else
            s(i,j)=b/a;
            end
        elseif  i==j+1
            s(i,j)=-1;
        elseif i==j-1
            s(i,j)=-1;
        else 
            s(i,j)=0;
        end
    end
end
...Create a coefficient matrix fot T_n
u=zeros(J-1);
for i=1:J-1
    for j=1:J-1
        if i==j
            if i==J-1
               u(i,j)=(c/a)+(1/(1+(h.*dx/k))); 
            else
            u(i,j)=c/a;
            end
        elseif  i==j+1
            u(i,j)=1;
        elseif i==j-1
            u(i,j)=1;
        else 
            u(i,j)=0;
        end
    end
end
%%Create matrices with known values
q=zeros(J-1,1);
for j=1:J-1
    if j==1
        q(j,1)=(d/a)+(2*T_0);
    elseif j==J-1
        q(j,1)=(d/a)+(2*h*dx*T_inf/(k+(h*dx)));
    else
       q(j,1)=d/a;
    end
end
...initial conditions
T_n=ones(J-1,1).*T_inf;
...Calculation of temperature distribution
for i=1:length(t)
T_n=linsolve(s,u*T_n+q);
end
T_end=(T_n(J-1,1)+(h*dx*T_inf/k))/(1+(h*dx/k));
T_n=[T_0;T_n;T_end]';
plot((0:J).*L./J,T_n,'o')
xlabel('Length (m)')
ylabel('Temperature (K)')