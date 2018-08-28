%% Input-State Linearization
clc
clear all

x0 = [.2 0];
tspan = [0 15];

[t,x] = ode45(@myode, tspan, x0);

xd = sin(t);
xd_d = cos(t);

figure(1);
plot(t,x(:,1)*180/pi,t,xd*180/pi)
legend("beta");

figure(3);
plot(t,x(:,2))
legend("betadot")


function dxdt = myode(t,x)

M = 1;
m = 0.1;
l = 0.5;
g = 9.81;
Fd = 1;

theta = x(1);
theta_dot = x(2);
xd = sin(t);
xd_d = cos(t);
xd_dd = -sin(t);

f = ((M+m)*g*sin(theta)-m*l*theta_dot^2*sin(theta)*cos(theta)+(M/m)*Fd*cos(theta))/((M+m)*l-m*l*cos(theta)^2);
b = cos(theta)/((M+m)*l-m*l*cos(theta)^2);

K = [137.78 25.97];

v = -K*[x(1)-xd; x(2)-xd_d];
u = -(f/g)+v;


dxdt = [theta_dot
        f + b*u];

        
end
