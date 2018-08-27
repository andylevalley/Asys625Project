%% Input-State Linearization
clc
clear all

x0 = [0 10*pi/180 0 0];
tspan = [0 15];

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


function dxdt = myode(t,x)

M = 1.0424;
m = 0.231;
l = 0.32;
It = 0.03155;
br = 0.00014;
beq = 9.582;
alpha = 0;
grav = 9.81;

X = [x(1); x(2); x(3); x(4)];
% Fc=0.5*tanh(100*x(3));
Fc = 0;
Msys = [M+m, -m*l*cos(x(2)-alpha); -m*l*cos(x(2)-alpha), It];
C = [beq, m*l*sin(x(2)-alpha)*x(4); 0, br];
D = [(M+m)*grav*sin(alpha); -m*grav*l*sin(x(2))]; 
f = Msys\(-C*[x(3);x(4)]-D - [Fc;0]);
g = Msys\([1;0]);

K = [-8.49 48.46 -15.38 9.67];
K = [0 48.46 0 9.67];

v = -K*X;

beta = x(2);
x_dot = x(3);
beta_dot = x(4);
u = v;
u = -f(2)/g(2) + v;

dxdt = [x(3)
        x(4)
        f(1) + g(1)*u
        f(2) + g(2)*u];

        
end
