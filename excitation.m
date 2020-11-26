clc 
clear
m=500;
n = 2;
a = 0.5;
b = 0.8;
p=10*eye(n,n);
e = wgn(m,1,0.1);
y = zeros(m,1);
phi = zeros(n,m);
theta=zeros(n,m);
% enable the desired input
% u = zeros(m,1);
% u = ones(m,1);
% u = wgn(m,1,1);
% u = 0.1*y(t);
% u =  0.1*sign(y(t));

for t=3:m
    y(t) = a * y(t-1) + b * u(t-1) + e(t);
    phi(:,t) = [y(t-1)  u(t-1)]';
    %RLS 
    p=p-p*(phi(:,t)*phi(:,t)')*p/(1+phi(:,t)'*p*phi(:,t));
    k = p *phi(:,t);
    theta(:,t)=theta(:,t-1)+k*(y(t)-phi(:,t)'*theta(:,t-1));
end

figure(1)
plot(1:m,y)
axis([0 m -6 6])
legend({'output: y'})
xlabel('time')
title('output: y')
figure(2)
plot(1:m,u)
hold on
plot(1:m,e)
axis([0 m -10 10])
legend({'Input: u','noise e'})
xlabel('time')
hold off 
figure (3) 
plot(1:m,theta(1,:))
hold on
plot(1:m,a*ones(1,m))
hold off
axis([0 m 0 3])
legend({'a estimated','a real'})
xlabel('time') 
figure (4)
plot(1:m,theta(2,:))
hold on
plot(1:m,b*ones(1,m))
hold off
axis([0 m 0 3])
legend({'b estimated','b real'})
xlabel('time')
 