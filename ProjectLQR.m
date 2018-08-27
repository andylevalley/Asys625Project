clc
clear all

x0 = [0 10*pi/180 0 0];
tspan = [0 30];

[t,x] = ode45(@myode, tspan, x0);

figure(1);
plot(t,x(:,1))
legend("position");

figure(2);
plot(t,x(:,3))
legend("velocity")

figure(3);
plot(t,x(:,2)*180/pi)
legend("beta")



function [system] = myode(t,x)

M = 1.0424;
m = 0.231;
l = 0.32;
It = 0.03155;
br = 0.00014;
beq = 9.582;
alpha = 10*pi/2;
g = 9.81;

gamma = (M+m)*It - m^2*l^2;
A = (1/gamma).*[0 0 gamma 0; 0 0 0 gamma; 0 m^2*l^2*g -It*beq -m*l*br; 0 (M+m)*m*g*l -m*l*beq -(M+m)*br];
B = (1/gamma).*[0;0;It;m*l];
X = [x(1); x(2); x(3); x(4)];

K = [-8.49, 48.46, -15.38, 9.67];

u = -K*X;


system = A*X + B*u;
dxdt = [x(3)
        x(4)
        system(3)
        system(4)];
end
                