%% Input-State Linearization

clc
clear all

x0 = [0 10*pi/180 0 0];
tspan = [0 30];

[t,x] = ode45(@myode_SMC, tspan, x0);

figure(1);
plot(t,x(:,1))
legend('eta');

figure(2);
plot(t,x(:,2)*180/pi)
legend('beta')






