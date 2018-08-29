%% Input-State Linearization
clc
clear all

x0 = [.2 0];
tspan = [0 15];
M = 1;
m = 0.1;
l = 0.5;
g = 9.81;
Fd = 10;

[t,x] = ode45(@myode, tspan, x0);

xd = sin(t);
xd_d = cos(t);
theta = x(:,1);
theta_dot = x(:,2);
f = (-(M+m).*g.*sin(theta)+m.*l.*theta_dot.^2.*sin(theta).*cos(theta)+(M/m).*Fd.*cos(theta))./(l.*M+m.*l.*cos(theta).^2);
b = cos(theta)./(l.*M+m.*l.*cos(theta).^2);

K = [137.78 25.97];
% 
% v = -K.*[x(:,1)-xd; x(:,2)-xd_d];
% u = -(f./b)+v;

figure(1);
plot(t,x(:,1)*180/pi,t,xd*180/pi)
legend("beta");

figure(3);
plot(t,x(:,2))
legend("betadot")


function dxdt = myode(t,x)

M = 1;
m = 0.2;
l = 0.32;
It = 0.03155;
br = 0.00014;
g = 9.81;
Fd = 10;

theta = x(1);
theta_dot = x(2);
xd = sin(t);
xd_d = cos(t);
xd_dd = pi+sin(t);

f = (-M*theta_dot*br - theta_dot*br*m + g*l*m^2*sin(theta) -...
    theta_dot^2*l^2*m^2*cos(theta)*sin(theta) ...
    + M*g*l*m*sin(theta) - (M/m)*Fd*cos(theta))/...
    (It*m + It*M - l^2*m^2*cos(theta)^2);
 
b = l*m*cos(theta)/(It*m + It*M - l^2*m^2*cos(theta)^2);

K = [138 27];

v = -K*[x(1)-xd; x(2)-xd_d];
u = -(f/b)+v;


dxdt = [theta_dot
        f + b*u];

        
end
