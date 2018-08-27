%% Input-State Linearization
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


function dxdt = myode(t,x)

M = 1.0424;
m = 0.231;
l = 0.32;
It = 0.03155;
br = 0.00014;
beq = 9.582;
alpha = 0;
alpha_hat = 5*pi/180;
eta = 0.1;
lambda1 = 2;
lambda2 = 2;
grav = 9.81;
Fc = 0;
% Fc=0.5*tanh(100*x(3));

X = [x(1); x(2); x(3); x(4)];
Msys = [M+m, -m*l*cos(x(2)-alpha); -m*l*cos(x(2)-alpha), It];
C = [beq, m*l*sin(x(2)-alpha)*x(4); 0, br];
D = [(M+m)*grav*sin(alpha); -m*grav*l*sin(x(2))]; 
f = Msys\(-C*[x(3);x(4)]-D - [Fc;0]);
g = Msys\([1;0]);

%% Calculate control law
M_hat = [M+m, -m*l*cos(x(2)-alpha_hat); -m*l*cos(x(2)-alpha_hat), It];
C_hat = [beq, m*l*sin(x(2)-alpha_hat)*x(4); 0, br];
D_hat = [(M+m)*grav*sin(alpha_hat); -m*grav*l*sin(x(2))];
f_hat = M_hat\(-C_hat*[x(3);x(4)]-D_hat - [Fc;0]);
g_hat = M_hat\([1;0]);

s1 = f(1) + 2*lambda1*x(3) + lambda1*x(1);
s2 = f(2) + 2*lambda2*x(4) + lambda2*x(2);
s_tot = s1+s2;
beta = 1;

u_hat = -f_hat(1) - f_hat(2) - lambda1*x(3) - lambda2*x(4);
F = (10*pi/180)*abs((f(1)+f(2)));
k = beta*(F + eta) + (beta - 1)*abs(u_hat);
u = (u_hat - k*sign(s_tot))/(g(1)+g(2));

dxdt = [x(3)
        x(4)
        f(1) + g(1)*u
        f(2) + g(2)*u];

        
end
