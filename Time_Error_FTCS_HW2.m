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
time=input('time= ');

...Number of pieces created on the length L
J=200;
dx=L/J;
W=10;
g=zeros(1,W-1);
G=zeros(1,W-1);
for w=0:W
    if w>1
       ...Midpoint temperature: 
       m=T(length(t),0.5*(J+2));
       g(1,w-1)=ta; 
    end
...Time step
ta=0.00001/w;
t=0:ta:time;
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
  if w>1
   G(1,w-1)=abs(T(length(t),0.5*(J+2))-m);
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