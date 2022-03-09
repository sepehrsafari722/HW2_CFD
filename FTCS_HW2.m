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
J=250;
dx=L/J;
...Time step
ta=0.00001;
t=0:ta:t;
...Calculate alfa and beta
alfa=k/(ro*C);
beta=(h*P)/(A*ro*C);

T=zeros(length(t),J+1);
T(1,:)=T_inf;
T(2:end,1)=T_0;

for i=2:length(t)
    for j=2:J
    T(i,j)=(alfa.*ta./dx.^2).*T(i-1,j-1)+((-2.*alfa.*ta./dx.^2)+1-(beta.*ta)).*T(i-1,j)+(alfa.*ta./dx.^2).*T(i-1,j+1)+(beta.*ta.*T_inf); 
       if j==J
          T(i,j+1)=(T(i,j)+(h*dx*T_inf/h))./(1+(h*dx/h));
       end
    end
end
plot(((0:J)./J).*L,T(length(t),:),'o')
xlabel('Length (m)')
ylabel('Temperature (K)')