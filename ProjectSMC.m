%% Sliding Mode Control
clc
clear all

M = 1;
m = 0.1;
l = 0.5;
m_hat = .01;
g = 9.81;
lambda = 10;
eta = 0.1;
phi = 0.1;


x0 = [0.2 0 eta/lambda];
tspan = [0 10];

[t,x] = ode45(@myode_SMC, tspan, x0);

xd = sin(t);
xd_d = cos(t);
xd_dd = -sin(t);

f = ((M+m).*g.*sin(x(:,1))-m.*l.*x(:,2).^2.*sin(x(:,1)).*cos(x(:,1)))./((M+m).*l-m.*l.*cos(x(:,1)).^2);
b = cos(x(:,1))./((4/3).*(M+m).*l-m.*l.*cos(x(:,1)).^2);


x_til = x(:,1) - xd;
x_til_dot = x(:,2) - xd_d;
% x_til_dotdot = f + xd_dd;

s = lambda*x_til + x_til_dot;

figure(1);
plot(t,x(:,1)*180/pi,t,xd*180/pi)
legend('beta');

figure(2);
plot(t,s,t,-x(:,3),t,x(:,3))
legend('s')

figure(3)
plot(t,x_til)
legend('xtil');






