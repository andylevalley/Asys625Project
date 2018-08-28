function XDOT = myode_Adaptive(t,X,param)
m = param.alpha1;
Fd = param.alpha2;
M = 1;
l = 0.5;
g = 9.81;


% f = ((M+m)*g*sin(theta)-m*l*theta_dot^2*sin(theta)*cos(theta)+Fd*cos(theta))/((M+m)*l-m*l*cos(theta)^2);
% b = cos(theta)/((M+m)*l-m*l*cos(theta)^2);


x = X(1);
xdot = X(2);
xdes = X(3);
xdesdot = X(4);
ahat1 = X(5);
ahat2 = X(6);
hhat = X(7);

lambda0 = param.lambda0;
gamma = param.gamma;
k = param.k;

r = -cos(t);



alpha1 = param.alpha1;
alpha2 = param.alpha2;


e = x-xdes;
edot = xdot-xdesdot;


sigma = xdot-xdesdot + lambda0*e;



f = ((M+m)*g*sin(x)-m*l*xdot^2*sin(x)*cos(x)+Fd*cos(x))/((M+m)*l-m*l*cos(x)^2);
f_alpha = ((M+alpha1)*g*sin(x)-alpha1*l*xdot^2*sin(x)*cos(x)+alpha2*cos(x))/((M+alpha1)*l-m*l*cos(x)^2);
b_alpha = cos(x)/((M+alpha1)*l-m*l*cos(x)^2);

xdesdotdot = r - f_alpha;
xrdotdot = xdesdotdot-lambda0*edot;

f_m = f - m*cos(x)/((M+m)*l-m*l*cos(x)^2);
f_Fd = (m*cos(x))/((M+m)*l-m*l*cos(x)^2);

ahat1dot = -gamma*sigma*f_m;
ahat2dot = -gamma*sigma*f_Fd;
hhatdot = -gamma*sigma*xrdotdot;

f_ahat = ((M+ahat1)*g*sin(x)-ahat1*l*xdot^2*sin(x)*cos(x)+ahat2*cos(x))/((M+ahat1)*l-m*l*cos(x)^2);


u = hhat*(xdesdotdot-lambda0*edot) -k*sigma + f_ahat;

xdotdot = b_alpha*u - f_alpha;

XDOT = [xdot; xdotdot; xdesdot; xdesdotdot; ahat1dot; ahat2dot; hhatdot];


end