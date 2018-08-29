clc
clear all

x0 = [0 0.2 0 0];
tspan = [0 15];

[t,x] = ode45(@myode, tspan, x0);

xd = sin(t);
xd_d = cos(t);

figure(1);
plot(t,x(:,2)*180/pi,t,xd*180/pi)
legend('beta');

figure(3);
plot(t,x(:,4))
legend('betadot')



function [system] = myode(t,x)

M = 1;
m = 0.2;
l = 0.32;
It = 0.03155;
br = 0.00014;
g = 9.81;
Fd = 1;

X = [x(1);x(2);x(3);x(4)];
gamma = (M+m)*It - m^2*l^2;
A = (1/gamma).*[0 1 gamma 0; 0 0 0 gamma; 0 m^2*l*g^2 0 -m*l*br; 0 (M+m)*m*g*l 0 (M+m)*br];
B = (1/gamma).*[0; 0; It; m*l];

K = [137.78 25.97];

xd = sin(t);
xd_d = cos(t);
x_til = x(2)-xd;
x_til_dot = x(4)-xd_d;

u = -K*[x_til; x_til_dot];


system = A*X + B*u + [0;0;0;1]*(M/m)*Fd*cos(x(2));
dxdt = system;
end
                