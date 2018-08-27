function dxdt = myode_SMC(t,x)

M = 1.0424;
m = 0.231;
l = 0.32;
It = 0.03155;
br = 0.00014;
beq = 9.582;
alpha = 0;
parm = .1;
alpha_hat = 0;
lambda = 10;
grav = 9.81;
Fc = 0;
% Fc=0.5*tanh(100*x(3));

% X = [x(1); x(2); x(3); x(4); x(5); x(6)];



beta = x(2);
beta_dot = x(4);
pos = x(1);
vel = x(3);
eta = pos - (It/(m*l))*log((1+sin(beta))/(cos(beta)));
eta_dot = vel - (It/(m*l*cos(beta)))*beta_dot;
eta_dotdot = -(grav+(It*beta_dot^2)/(m*l*cos(beta)))*tan(beta)+(br*beta_dot)/(m*l*cos(beta));


Msys = [M+m, -m*l*cos(x(2)-alpha); -m*l*cos(x(2)-alpha), It];
C = [beq, m*l*sin(x(2)-alpha)*x(4); 0, br];
D = [(M+m)*grav*sin(alpha); -m*grav*l*sin(x(2))]; 
f = Msys\(-C*[x(3);x(4)]-D-[Fc;0]);
g = Msys\[1;0];

s = eta_dot + lambda*eta;
k = 0.1;
ueq = -eta_dotdot - lambda*eta_dot;
u = ueq-k*sign(s);


dxdt = [vel
        beta_dot
        f(1) + g(1)*u
        f(2) + g(2)*u];
    
end