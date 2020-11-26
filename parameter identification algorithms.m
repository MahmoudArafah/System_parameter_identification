%Recursive Least Square
clc 
clear
m=1000;
n = 2;
a = 0.7;
c = -0.5;
p=10*eye(n,n);
e = wgn(m,1,.1);
y = zeros(m,1);
phi = zeros(n,m);
theta=zeros(n,m);
for t=3:m
    y(t) = a * y(t-1) + c * e(t-1) + e(t) ;
    phi(:,t) = [y(t-1) 0]';
    p=p-p*(phi(:,t)*phi(:,t)')*p/(1+phi(:,t)'*p*phi(:,t));
    k = p *phi(:,t);
    theta(:,t)=theta(:,t-1)+k*(y(t)-phi(:,t)'*theta(:,t-1));
end
figure(1)
plot(1:m,y)
axis([0 m -6 6])
legend({'output: y'})
xlabel('time')
title('output: y using RLS')
figure(2)
plot(1:m,e)
axis([0 m -10 10])
legend({'e'})
xlabel('time')
title('e using RLS')
figure (3) 
plot(1:m,theta(1,:))
hold on
plot(1:m,a*ones(1,m))
hold off
axis([0 m -1.5 1.5])
legend({'a estimated','a real'})
xlabel('time') 
title('estimation of a using RLS')
figure (4)
plot(1:m,theta(2,:))
hold on
plot(1:m,c*ones(1,m))
hold off
legend({'c estimated','c real'})
xlabel('time')
title('estimation of c using RLS')


% Extended Least Square
clc 
clear
m=1000;
n = 3;
a = 0.7;
c = -0.5;
p=10*eye(n,n);
e = wgn(m,1,1);
y = zeros(m,1);
phi = zeros(n,m);
theta=zeros(n,m);
epsilon = zeros(m,1);
for t=3:m
    y(t) = a * y(t-1) + c * e(t-1) + e(t);
    phi(:,t) = [y(t-1)  e(t-1) epsilon(t-1)]';
    epsilon(t) = y(t) - phi(:,t)' * theta(:,t-1);
    p=p-p*(phi(:,t)*phi(:,t)')*p/(1+phi(:,t)'*p*phi(:,t));
    theta(:,t)=theta(:,t-1)+p*phi(:,t)*(y(t)-phi(:,t)'*theta(:,t-1));
end
figure(5)
plot(1:m,y)
axis([0 m -6 6])
legend({'output: y'})
xlabel('time')
title('output: y using ELS')
figure(6)
plot(1:m,e)
axis([0 m -10 10])
legend({'e'})
xlabel('time')
title('e using ELS')
figure (7) 
plot(1:m,theta(1,:))
hold on
plot(1:m,a*ones(1,m))
hold off
legend({'a estimated','a real'})
xlabel('time') 
title('estimation of a using ELS')
figure (8)
plot(1:m,theta(2,:))
hold on
plot(1:m,c*ones(1,m))
hold off
legend({'c estimated','c real'})
xlabel('time')
title('estimation of c using ELS')

% Approximate Maximum Likelihood
clc 
clear
m=1000;
n = 3;
a = 0.7;
c = -0.5;
p=10*eye(n,n);
e = wgn(m,1,1);
y = zeros(m,1);
phi = zeros(n,m);
theta=zeros(n,m);
eta = zeros(m,1);
for t=3:m
   y(t) = a * y(t-1) + c * e(t-1) + e(t);
   phi(:,t) = [y(t-1)  e(t-1) eta(t-1)]';
   eta(t) = y(t) - phi(:,t)' * theta(:,t);
   p=p-p*(phi(:,t)*phi(:,t)')*p/(1+phi(:,t)'*p*phi(:,t));
   theta(:,t)=theta(:,t-1)+p*phi(:,t)*(y(t)-phi(:,t)'*theta(:,t-1));
end
figure(9)
plot(1:m,y)
axis([0 m -6 6])
legend({'output: y'})
xlabel('time')
title('output: y using AML')
figure(10)
plot(1:m,e)
axis([0 m -10 10])
legend({'e'})
xlabel('time')
title('e using AML')
figure (11) 
plot(1:m,theta(1,:))
hold on
plot(1:m,a*ones(1,m))
hold off
legend({'a estimated','a real'})
xlabel('time') 
title('estimation of a using AML')
figure (12)
plot(1:m,theta(2,:))
hold on
plot(1:m,c*ones(1,m))
hold off
legend({'c estimated','c real'})
xlabel('time')
title('estimation of c using AML')