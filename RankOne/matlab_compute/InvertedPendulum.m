clc
clear
syms theta(t) x(t)
syms N P M m I_m L g F T
syms theta_ddot x_ddot theta_dot x_dot

% M = 1.87500000;
% g = 9.805;
% 
% m = 1.0000;
% I_m = 0.01354167;
% L = 0.4;
% l = 0.2;

l = L / 2;

funcs = [
    M * diff(x, t, t) == F - N;
    m * diff(-l * cos(theta), t, t) == P - m * g;
    N == diff(m * (l * sin(theta) + x), t, t);
    I_m * diff(theta, t, t) == P * l * sin(theta) - N * l * cos(theta);
    ];

funcs = subs(funcs, [diff(diff(theta, t), t), diff(diff(x, t), t), diff(theta, t), diff(x, t)], [theta_ddot, x_ddot, theta_dot, x_dot]);

[theta_ddot, x_ddot, P, N] = solve(funcs, [theta_ddot, x_ddot, P, N]);
X_dot = [theta_dot, theta_ddot, x_dot, x_ddot];
X = [theta(t), theta_dot, x(t), x_dot];
U = [F];


A = jacobian(X_dot, X);
B = jacobian(X_dot, U);

A = subs(A, [theta(t), theta_dot, x(t), x_ddot, F], [0, 0, 0, 0, 0]);
B = subs(B, [theta(t), theta_dot, x(t), x_ddot, F], [0, 0, 0, 0, 0]);

A = double(A);
B = double(B);

C = eye(4);
D = zeros(4,1);
Q = diag([100 1 100 1]);
R = diag([1]);
sys = ss(A, B, C, D);
KLQR = lqr(sys, Q, R);%得到反馈增益矩阵
