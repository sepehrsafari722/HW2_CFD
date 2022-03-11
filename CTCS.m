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
J=35;
dx=L/J;
...Time step
ta=0.009;
t=0:ta:t;
...Calculate alfa and beta
alfa=k/(ro*C);
beta=(h*P)/(A*ro*C);
...Calculation of coefficients a, b, d
a=-(alfa.*ta./dx.^2);
b=((2.*alfa.*ta./dx.^2)+1+(beta.*ta));
d=-(beta.*ta.*T_inf);
...Create a coefficient matrix
s=zeros(J-1);
for i=1:J-1
    for j=1:J-1
        if i==j
            if i==J-1
               s(i,j)=b+(a/(1+(h.*dx/k))); 
            else
            s(i,j)=b;
            end
        elseif  i==j+1
            s(i,j)=a;
        elseif i==j-1
            s(i,j)=a;
        else 
            s(i,j)=0;
        end
    end
end
%%Create matrices with known values
q=zeros(J-1,1);
for j=1:J-1
    if j==1
        q(j,1)=-d-(a*T_0);
    elseif j==J-1
        q(j,1)=-a*(h*dx*T_inf/(k+(h*dx)))-d;
    else
       q(j,1)=-d;
    end
end
...initial conditions
T_n=ones(J-1,1).*T_inf;
...Calculation of temperature distribution
T_nplus1=linsolve(s,T_n+q);
T_end=(T_nplus1(J-1,1)+(h*dx*T_inf/k))/(1+(h*dx/k));
T_nplus1=[T_0;T_nplus1;T_end]';
z=zeros(length(t),J+1);
z(1,:)=T_inf;
z(2,:)=T_nplus1;
z(3:end,1)=T_0;
for  i=3:length(t)
for  j=2:J 
    z(i,j)=(2*alfa*ta/(dx^2))*z(i-1,j-1)-((4*alfa*ta/(dx^2))+(2*beta*ta))*z(i-1,j)+(2*alfa*ta/(dx^2))*z(i-1,j+1)+z(i-2,j)+2*(beta*ta*T_inf); 
    if j==J
       z(i,j+1)=(z(i,j)+(h*dx*T_inf/k))/(1+(h*dx/k));
    end
end
end
plot((0:J).*L./J,z(length(t),:),'o')
xlabel('Length (m)')
ylabel('Temperature (K)')