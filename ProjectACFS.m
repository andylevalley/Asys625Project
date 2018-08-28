clc; clear all; close all;
%adaptive control, nonlinear system, example (mass-spring-damper)
%ASYS 625 SU18, Hess (use at own risk)
x0 = 0; xdot0 = 0; %plant ICs
xdes0 = 0; xdesdot0 = 0; %xdes ICs
t0 = 0; tf = 20; tspan = [t0 tf]; %time span
%% desired params
m = .1;
Fd = 1;
param.alpha1 = m; 
param.alpha2 = Fd;
%%
%assumes h = 1
ahat10 = 0.5;
ahat20 = .5;
hhat0 = 1;

%% tuning
param.gamma = 1;
param.lambda0 = 3;
param.k = 2;

%%

X0 = [x0; xdot0; xdes0; xdesdot0; ahat10; ahat20; hhat0];
options = [];
[t,X] = ode45(@myode_Adaptive, tspan, X0, options, param);


x = X(:,1); xdot = X(:,2); xdes = X(:,3); xdesdot = X(:,4);
ahat1 = X(:,5); ahat2 = X(:,6); hhat = X(:,7);

%% Plot

r = sin(t);
alpha1 = param.alpha1;
alpha2 = param.alpha2;

figure;
subplot(221);
plot(t,x,t,xdes); legend('$x$','$x_d$')
xlabel('$t(s)$'); ylabel('State and Desired');
subplot(222);
plot(t,alpha1*ones(size(t)), 'b', t, ahat1,'b-.',...
     t,alpha2*ones(size(t)), 'r', t, ahat2,'r-.');
legend('$\alpha_1$','$\hat{a}_1$',...
       '$\alpha_2$','$\hat{a}_2$');
xlabel('$t(s)$'); ylabel('Parameters');
subplot(223);
plot(t,x-xdes); 
xlabel('$t(s)$'); ylabel('Error ($x-x_d$)');
