clc
clear
syms theta(t) x(t) phi(t)
syms N1 P1 m1 L1 I1 N2 P2 m2 L2 I2 M g F
syms theta_ddot x_ddot theta_dot x_dot phi_ddot phi_dot

M = 1.87500000;
g = 9.805;

m1 = 1.0000;
I1 = 0.01354167;
L1 = 0.4;

m2 = 1.0000;
I2 = 0.01354167;
L2 = 0.4;

funcs = [
    F - N1 == M * diff(diff(x, t), t);
    N1 - N2 == m1 * diff(diff(x + L1 * sin(theta) / 2, t), t);
    P1 - P2 - m1 * g == m1 * diff(diff(L1 * cos(theta) / 2, t), t);
    I1 * diff(diff(theta, t), t) == (P1 + P2) * sin(theta) * L1 / 2 - (N1 + N2) * cos(theta) * L1 / 2;
    N2 == m2 * diff(diff(x + L1 * sin(theta) / 2 + L2 * sin(phi) / 2, t), t);
    P2 - m2 * g == m2 * diff(diff(L1 * cos(theta) + L2 * cos(phi) / 2, t), t);
    I2 * diff(diff(phi, t), t) == P2 * L2 * sin(phi) / 2 - N2 * L2 * cos(phi) / 2
    ];


funcs = subs(funcs, ...
    [diff(diff(theta, t), t), diff(diff(x, t), t), diff(diff(phi, t), t), diff(theta, t), diff(x, t), diff(phi, t)], ...
    [theta_ddot, x_ddot, phi_ddot, theta_dot, x_dot, phi_dot]);

[theta_ddot, x_ddot, phi_ddot, N1, P1, N2, P2] = solve(funcs, [theta_ddot, x_ddot, phi_ddot, N1, P1, N2, P2]);
theta_ddot = simplify(theta_ddot);
x_ddot = simplify(x_ddot);
phi_ddot = simplify(phi_ddot);

X_dot = [theta_dot, theta_ddot, x_dot, x_ddot, phi_dot, phi_ddot];
X = [theta(t), theta_dot, x(t), x_dot, phi(t), phi_dot];
U = [F];

A = jacobian(X_dot, X);
B = jacobian(X_dot, U);

A = subs(A, [theta(t), theta_dot, x_ddot, phi(t), phi_dot, F], [0, 0, 0, 0, 0, 0]);
B = subs(B, [theta(t), theta_dot, x_ddot, phi(t), phi_dot, F], [0, 0, 0, 0, 0, 0]);

A = double(A);
B = double(B);

C = eye(6);
D = zeros(6, 1);
Q = diag([100 50 100 10 5000 50]);
R = diag([1]);
sys = ss(A, B, C, D);
KLQR = lqr(sys, Q, R)%得到反馈增益矩阵

