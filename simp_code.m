clc
clear
close all
%Numerical solution of  Volterra Integral second Kind equation
%by simpson method
% Eq=> y(x)=f(x)+int(k(x,t)*y(t)dy,a,x)
%% input parameters
a=input('inter lower limit of int a=');
b=input('inter upper limit of int instead of x=');
   f_x=input('inter known function in form  @(x) f(x)        :');
   k_xt=input('inter equation kernal in form @(x,t) k(x,t)   :');
             n=20;h=(b-a)/n;                %the number of interval divisions that must be even.
%% solution
x=linspace(a,b,n+1);
fvec=f_x(x);
uvec=zeros(1,n+1);
uvec(1)=fvec(1);
for i=1:n
    uvec(i+1)=uvec(i);
    kvec=k_xt(x(i+1),x(1:i+1)).*uvec(1:i+1);
    uvec(i+1)=fvec(i+1)+(h/3)*((2*sum(kvec(2:2:i))+4*sum(kvec(3:2:i)))+(kvec(1)+kvec(i+1)));  
end
u=uvec;
%% comparsion with the exact suolution
x=linspace(a,b,n+1);
u_ex=input('inter exact solution of equation like (x) x.*(x-sin(x))./(x.^2):');
y=abs(u_ex-u);
m=[x',u',u_ex',y'];
writematrix(m,'Data2.xlsx','Range','A2:D22')
plot(x,u,'*',x,u_ex,'r')
grid on
legend('Approximation Solution','Exact Solution')
xlabel('x')
ylabel('u(x)')
title('Simpson method')
