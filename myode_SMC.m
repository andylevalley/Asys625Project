function dxdt = myode_SMC(t,x)

M = 1;
m = 0.1;
l = 0.5;
m_hat = .1;
g = 9.81;
lambda = 10;
eta = 0.1;
Fd = 1;
Fd_hat = 1;

theta = x(1);
theta_dot = x(2);
phi = x(3);

f = ((M+m)*g*sin(theta)-m*l*theta_dot^2*sin(theta)*cos(theta)+(M/m)*Fd*cos(theta))/((M+m)*l-m*l*cos(theta)^2);
b = cos(theta)/((M+m)*l-m*l*cos(theta)^2);
f_hat = ((M+m_hat)*g*sin(theta)-m_hat*l*theta_dot^2*sin(theta)*cos(theta)+(M/m)*Fd_hat*cos(theta))/((M+m_hat)*l-m_hat*l*cos(theta)^2);
b_hat = cos(theta)/((M+m_hat)*l-m*l*cos(theta)^2);


xd = sin(t);
xd_d = cos(t);
xd_dd = -sin(t);

x_til = x(1) - xd;
x_til_dot = x(2) - xd_d;
x_til_dotdot = f + xd_dd;

s = lambda*x_til + x_til_dot;

u_hat = -f_hat + xd_dd - lambda*x_til_dot;
F = 3*abs(f);
k = F+eta;

f_d = ((M+m)*g*sin(xd)-m*l*xd_d^2*sin(xd)*cos(xd)+(M/m)*Fd*cos(xd))/((M+m)*l-m*l*cos(xd)^2);
k_d = 3*abs(f_d);

k_bar = k - k_d + lambda*phi;
phi_dot = -lambda*phi + k_d;

if abs(s) <= phi
    u = b_hat^(-1)*(u_hat - k_bar*(s/phi));
else
    u = b_hat^(-1)*(u_hat - k_bar*sign(s));
end

dxdt = [theta_dot
        f + b*u
        phi_dot];
    
end

