%% Input-State Linearization
clc
clear all

x0 = [.2 0];
tspan = [0 15];
M = 2.4;
m = 0.23;
l = 0.36;
g = 9.81;
Fd = 1;

[t,x] = ode45(@myode, tspan, x0);

xd = sin(t);
xd_d = cos(t);
theta = x(:,1);
theta_dot = x(:,2);
f = (-(M+m).*g.*sin(theta)+m.*l.*theta_dot.^2.*sin(theta).*cos(theta)+(M/m).*Fd.*cos(theta))./(l.*M+m.*l.*cos(theta).^2);
b = cos(theta)./(l.*M+m.*l.*cos(theta).^2);

K = [40 8];

v = -40*(x(:,1)-xd) - 8*(x(:,2)-xd_d);
u = -(f./b)+v;

figure(1);
plot(t,x(:,1)*180/pi,t,xd*180/pi)
legend("beta");

figure(3);
plot(t,x(:,2))
legend("betadot")


function dxdt = myode(t,x)

M = 2.4;
m = 0.23;
l = 0.36;
g = 9.81;
g = 9.81;
Fd = 10;

theta = x(1);
theta_dot = x(2);
xd = sin(t);
xd_d = cos(t);
xd_dd = -sin(t);

f = (-(M+m)*g*sin(theta) + m*l*cos(theta)*sin(theta)*theta_dot^2-(M/m)*Fd*cos(theta))/(m*l*cos(theta)^2-(M+m)*l);
b = cos(theta)/(m*l*cos(theta)^2-(M+m)*l);

K = [-137.78 -25.97];

v = -K*[x(1)-xd; x(2)-xd_d];
u = -(f/b)+v;
% u = v;


dxdt = [theta_dot
        f + b*u];

        
end
