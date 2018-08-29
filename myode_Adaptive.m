function XDOT = myode_Adaptive(t,X,param)

M = 1;
m = 0.2;
l = 0.32;
It = 0.03155;
br = 0.00014;
g = 9.81;
Fd = param.alpha1;


x = X(1);
xdot = X(2);
xdes = X(3);
xdesdot = X(4);
ahat1 = X(5);
hhat = X(6);

lambda0 = param.lambda0;
gamma = param.gamma;
k = param.k;


xdes = sin(t);
xdesdot = cos(t);



alpha1 = param.alpha1;

xrefdot = xdesdot - lambda*
s = xdot-xrefdot;


f_alpha = (-M*xdot*br - xdot*br*m + g*l*m^2*sin(x) -...
    xdot^2*l^2*m^2*cos(x)*sin(x) ...
    + M*g*l*m*sin(x) - (M/m)*alpha1*cos(x))/...
    (It*m + It*M - l^2*m^2*cos(x)^2);

b_alpha = l*m*cos(x)/(It*m + It*M - l^2*m^2*cos(x)^2);
xdesdotdot = r - f_alpha;

xrdotdot = xdesdotdot-lambda0*edot;

f_Fd = -((M/m)*cos(x))/(It*m + It*M - l^2*m^2*cos(x)^2);
ahat1dot = -gamma*sigma*f_Fd;
hhatdot = -gamma*sigma*xrdotdot;

f_ahat = (-M*xdot*br - xdot*br*m + g*l*m^2*sin(x) -...
    xdot^2*l^2*m^2*cos(x)*sin(x) ...
    + M*g*l*m*sin(x) - (M/m)*ahat1*cos(x))/...
    (It*m + It*M - l^2*m^2*cos(x)^2);

u = hhat*(xdesdotdot-lambda0*edot) -k*sigma + f_ahat;


xdotdot = f_alpha + b_alpha*u;

XDOT = [xdot; xdotdot; xdesdot; xdesdotdot; ahat1dot; hhatdot];

end