clc; clear all; close all;
%adaptive control, nonlinear system, example (mass-spring-damper)

x0 = 0.2; xdot0 = 0; %plant ICs
xdes0 = 0; xdesdot0 = 0; %xdes ICs
t0 = 0; tf = 30; tspan = [t0 tf]; %time span
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
param.gamma = 20;
param.lambda0 = 100;
param.k = 2;

%%

X0 = [x0; xdot0; xdes0; xdesdot0; ahat10; hhat0];
options = [];
[t,X] = ode45(@myode_MRAC, tspan, X0, options, param);


x = X(:,1); xdot = X(:,2); xdes = X(:,3); xdesdot = X(:,4);
ahat1 = X(:,5); hhat = X(:,6);

%% Plot

r = sin(t);
alpha1 = param.alpha1;

figure(1)
plot(t,x*180/pi,t,xdes*180/pi); legend('$x$','$x_d$')
xlabel('$t(s)$'); ylabel('State and Desired');

figure(2)
plot(t,alpha1*ones(size(t)), 'b', t, ahat1,'b-.');
legend('$\alpha_1$','$\hat{a}_1$');
xlabel('$t(s)$'); ylabel('Parameters');

figure(3)
plot(t,x-xdes); 
xlabel('$t(s)$'); ylabel('Error ($x-x_d$)');

function XDOT = myode_MRAC(t,X,param)
M = 2.4;
m = 0.23;
l = 0.36;
g = 9.81;
zeta1 = 27.97;
zeta2 = 137.78;

x = X(1);
xdot = X(2);
x_model = X(3);
xdot_model = X(4);
ahat1 = X(5);
hhat = X(6);

lambda0 = param.lambda0;
gamma = param.gamma;
k = param.k;

r = sin(t);

xdotdot_model = -zeta1*xdot_model - zeta2*x_model - zeta2*r;

e = x - x_model;
edot = xdot-xdot_model;
xdot_ref = xdot_model - lambda0*e;
xdotdot_ref = xdotdot_model - lambda0*edot;

s = xdot - xdot_ref;

% f = (-(M+m)*g*sin(theta) + m*l*cos(theta)*sin(theta)*theta_dot^2-(M/m)*Fd*cos(theta))/(m*l*cos(theta)^2-(M+m)*l);
% b = cos(theta)/(m*l*cos(theta)^2-(M+m)*l);
% f_alpha1 = (-(M+m)*g*sin(x) + m*l*cos(x)*sin(x)*xdot^2-(M/m)*alpha1*cos(x))/(m*l*cos(x)^2-(M+m)*l);
% b_alpha1 = cos(x)/(m*l*cos(x)^2-(M+m)*l);

f_ahat1 = (-(M+m)*g*sin(x) + m*l*cos(x)*sin(x)*xdot^2-(M/m)*ahat1*cos(x))/(m*l*cos(x)^2-(M+m)*l);
b_ahat1 = cos(x)/(m*l*cos(x)^2-(M+m)*l);
f_Fb = (-(M/m)*cos(x))/(m*l*cos(x)^2-(M+m)*l);

u = (hhat/b_ahat1)*xdotdot_ref -k*s + (ahat1/b_ahat1)*f_Fb;
ahat1dot = -gamma*s*f_Fb;
hhatdot = -gamma*s*xdotdot_ref;

xdotdot = f_ahat1 + b_ahat1*u;

XDOT = [xdot; xdotdot; xdot_model; xdotdot_model; ahat1dot; hhatdot];


end