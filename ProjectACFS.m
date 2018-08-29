clc; clear all; close all;
%adaptive control, nonlinear system, example (mass-spring-damper)

x0 = pi; xdot0 = 0; %plant ICs
xdes0 = 0; xdesdot0 = 0; %xdes ICs
t0 = 0; tf = 100; tspan = [t0 tf]; %time span
%% desired params
m = .1;
Fd = 1;
% param.alpha1 = m; 
param.alpha1 = Fd;
%%
%assumes h = 1
ahat10 = 0.5;
hhat0 = 1;

%% tuning
param.gamma = 2;
param.lambda0 = 500;
param.k = 2;

%%

X0 = [x0; xdot0; xdes0; xdesdot0; ahat10; hhat0];
options = [];
[t,X] = ode45(@myode_Adaptive, tspan, X0, options, param);


x = X(:,1); xdot = X(:,2); xdes = X(:,3); xdesdot = X(:,4);
ahat1 = X(:,5); hhat = X(:,6);

%% Plot

r = sin(t);
alpha1 = param.alpha1;

figure;
subplot(221);
plot(t,x*180/pi,t,xdes*180/pi); legend('$x$','$x_d$')
xlabel('$t(s)$'); ylabel('State and Desired');
subplot(222);
plot(t,alpha1*ones(size(t)), 'b', t, ahat1,'b-.');
legend('$\alpha_1$','$\hat{a}_1$');
xlabel('$t(s)$'); ylabel('Parameters');
subplot(223);
plot(t,x-xdes); 
xlabel('$t(s)$'); ylabel('Error ($x-x_d$)');
