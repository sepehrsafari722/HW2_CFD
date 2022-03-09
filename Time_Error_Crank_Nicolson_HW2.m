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

W=10;
g=zeros(1,W-1);
G=zeros(1,W-1);
for w=0:W
    if w>1
       ...Midpoint temperature: 
       m=T_n(0.5*(J+2),1);
       g(1,w-1)=ta; 
    end
...Number of pieces created on the length L
J=300;
dx=(L/J);
...Time step
ta=0.1+(w*10^-2);
t=0:ta:t;
...Calculate alfa and beta
alfa=k/(ro*C);
beta=(h*P)/(A*ro*C);
...Calculation of coefficients a, b, d
a=(alfa*ta/(2*dx^2));
b=((alfa*ta/dx^2))+1;
c=((-alfa*ta/dx^2)+1-(beta*ta));
d=beta*ta*T_inf;
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
...Calculation of successive error
if w>1
G(1,w-1)=abs(T_n(0.5*(J+2),1)-m);
end
end
figure(1)
loglog(g,G,'o')
xlabel('ta (s)')
ylabel('Error (K)')
figure(2)
slope=gradient(log(G))./gradient(log(g));
plot(g,slope,'o')
xlabel('ta (s)')
ylabel('Slope')