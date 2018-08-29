clc
clear all

x0 = [0.2 0 ];
tspan = [0 15];

[t,x] = ode45(@myode, tspan, x0);

xd = sin(t);
xd_d = cos(t);

figure(1);
plot(t,x(:,1)*180/pi,t,xd*180/pi)
legend('beta');

figure(3);
plot(t,x(:,2))
legend('betadot')



function [dxdt] = myode(t,x)

M = 2.4;
m = 0.23;
l = 0.36;
g = 9.81;
Fd = 1;

X = [x(1);x(2)];
A = [0 1; (M+m)*g/(M*l) 0];
B = [0; -1/(M*l)];
D = [0; -1/(m*l)];


K = [-137.78 -25.97];

xd = sin(t);
xd_d = cos(t);
x_til = x(1)-xd;
x_til_dot = x(2)-xd_d;

u = -K*[x_til; x_til_dot];


system = A*X + B*u + D*Fd;
dxdt = system;
end
                