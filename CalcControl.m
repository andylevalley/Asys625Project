function u = CalcControl(t,X,param)
gamma = param.gamma;
P = param.P;
alpha2 = param.alpha2;
alpha1 = param.alpha1;
alpha0 = param.alpha0;
beta1 = param.beta1;
beta0 = param.beta0;
%%
x = X(:,1);
xdot = X(:,2);
ahat0 = X(:,3);
ahat1 = X(:,4);
ahat2 = X(:,5);
xm = X(:,6);
xmdot = X(:,7);
%%
if param.ref == 0
    r = 1;
elseif param.ref == 1
    r = 4*sin(t);
elseif param.ref == 2
    r = 4*sin(t) + 3*cos(t);
elseif param.ref == 3
   r = 4*sin(t) + 3*cos(t) + 4*sin(2*t);
end
edot = xdot-xmdot;
e = x-xm;

xmdotdot = (1/alpha2)*(r - alpha1*xmdot-alpha0*xm);
z = xmdotdot - beta1*edot - beta0*e;
v = [z; xdot; x];

u = z.*ahat2 + xdot.*ahat1 + x.*ahat0;

b = [0;1];
evec = [e; edot];



end