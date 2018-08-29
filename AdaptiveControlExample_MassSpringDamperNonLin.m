clc; clear all; close all;
%adaptive control, nonlinear system, example (mass-spring-damper)
%ASYS 625 SU18, Hess (use at own risk)
x0 = 0; xdot0 = 0; %plant ICs
xdes0 = 0; xdesdot0 = 0; %xdes ICs
t0 = 0; tf = 20; tspan = [t0 tf]; %time span
%% desired params
b = 0.7; %damping
k0 = 1; %linear spring term
k1 = 0.2; %cubic spring term
m = 1;  %mass
alpha1 = (b/m); 
alpha2 = (k0/m);
alpha3 = (k1/m);
%%
%assumes h = 1
ahat10 = 1;
ahat20 = 1;
ahat30 = 1;
hhat0 = 1;
%% tuning
param.gamma = 2;
param.lambda0 = 2;
param.k = 2;
%%
param.ref = 2;
%1 for step, 2 for sin(t), 3 for sin(t) + sin(2t), 4 for
%sin(t)+sin(2t)+sin(4t)
param.const_ref = 1;
param.b = b;
param.m = m;
param.k0 = k0;
param.k1 = k1;

X0 = [x0;xdot0; xdes0; xdesdot0; ahat10; ahat20; ahat30; hhat0];
options = [];
[t,X] = ode45(@ExampleEOMs, tspan, X0, options, param);

r = CalcRef(t,param);

x = X(:,1); xdot = X(:,2); xdes = X(:,3); xdesdot = X(:,4);
ahat1 = X(:,5); ahat2 = X(:,6); ahat3 = X(:,7); hhat = X(:,8);

figure;
subplot(221);
plot(t,x,t,xdes); legend('$x$','$x_d$')
xlabel('$t(s)$'); ylabel('State and Desired');
subplot(222);
plot(t,alpha1*ones(size(t)), 'b', t, ahat1,'b-.',...
     t,alpha2*ones(size(t)), 'r', t, ahat2,'r-.',...
     t,alpha3*ones(size(t)), 'k', t, ahat3,'k-.');
legend('$\alpha_1$','$\hat{a}_1$',...
       '$\alpha_2$','$\hat{a}_2$',...
       '$\alpha_3$','$\hat{a}_3$');
xlabel('$t(s)$'); ylabel('Parameters');
subplot(223);
plot(t,x-xdes); 
xlabel('$t(s)$'); ylabel('Error ($x-x_d$)');
   
[u,sigma] = CalcControl(t,X,param);
subplot(224);
plot(t,u,t,sigma); 
xlabel('$t(s)$'); ylabel('$u$ and $\sigma$');
legend('$u$','$\sigma$')

figure;
plot(x-xdes,xdot-xdesdot,'-.',x(1)-xdes(1),xdot(1)-xdesdot(1),'o',x(end)-xdes(end),xdot(end)-xdesdot(end),'s');
xlabel('$e$'); ylabel('$\dot{e}$');

%%    
function XDOT = ExampleEOMs(t,X,param)
x = X(1);
xdot = X(2);
xdes = X(3);
xdesdot = X(4);
ahat1 = X(5);
ahat2 = X(6);
ahat3 = X(7);
hhat = X(8);

lambda0 = param.lambda0;
gamma = param.gamma;
k = param.k;

b = param.b;
k0 = param.k0;
k1 = param.k1;
m = param.m;

r = CalcRef(t,param);



alpha1 = (b/m);
alpha2 = (k0/m);
alpha3 = (k1/m);

e = x-xdes;
edot = xdot-xdesdot;


sigma = xdot-xdesdot + lambda0*e;



f1 = xdot*abs(xdot);
f2 = x;
f3 = x^3;
xdesdotdot = r - alpha1*f1 - alpha2*f2 - alpha3*f3;

xrdotdot = xdesdotdot-lambda0*edot;

ahat1dot = -gamma*sigma*f1;
ahat2dot = -gamma*sigma*f2;
ahat3dot = -gamma*sigma*f3;
hhatdot = -gamma*sigma*xrdotdot;


u = hhat*(xdesdotdot-lambda0*edot) -k*sigma + ahat1*f1 + ahat2*f2 + ahat3*f3;

xdotdot = u - alpha1*f1 - alpha2*f2 - alpha3*f3;

XDOT = [xdot; xdotdot; xdesdot; xdesdotdot; ahat1dot; ahat2dot; ahat3dot; hhatdot];



end


%%
function [u,sigma] = CalcControl(t,X,param)
x = X(:,1);
xdot = X(:,2);
xdes = X(:,3);
xdesdot = X(:,4);
ahat1 = X(:,5);
ahat2 = X(:,6);
ahat3 = X(:,7);
hhat = X(:,8);

lambda0 = param.lambda0;
gamma = param.gamma;
k = param.k;

b = param.b;
k0 = param.k0;
k1 = param.k1;
m = param.m;


r = CalcRef(t,param);

alpha1 = (b/m);
alpha2 = (k0/m);
alpha3 = (k1/m);

e = x-xdes;
edot = xdot-xdesdot;


sigma = xdot-xdesdot + lambda0*e;



f1 = xdot.*abs(xdot);
f2 = x;
f3 = x.^3;
xdesdotdot = r - alpha1*f1 - alpha2*f2 - alpha3*f3;
xrdotdot = xdesdotdot-lambda0*edot;

ahat1dot = -gamma*sigma.*f1;
ahat2dot = -gamma*sigma.*f2;
ahat3dot = -gamma*sigma.*f3;
hhatdot = -gamma*sigma.*xrdotdot;


u = hhat.*(xdesdotdot-lambda0.*edot) -k.*sigma + ahat1.*f1 + ahat2.*f2 + ahat3.*f3;

xdotdot = u - alpha1*f1 - alpha2*f2 - alpha3*f3;

XDOT = [xdot; xdotdot; xdesdot; xdesdotdot; ahat1dot; ahat2dot; ahat3dot; hhatdot];



end


function r = CalcRef(t,param)
if param.ref == 1
    r = param.const_ref;
elseif param.ref == 2
    r = sin(t);
elseif param.ref == 3
    r = sin(t) + sin(2*t);
elseif param.ref == 4
    r = sin(t) + sin(2*t) + sin(4*t);
end
end



