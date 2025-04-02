predict_n = 10;
predict_step = 0.04;

x0 = [0; 0; 0; 0];
x_ref = [0; 0; 1; 0];
F_val = 0;
X = x0;
J = 0;
for i = 1 : 1 : predict_n
    theta_ddot_val = subs(theta_ddot, [theta(t), theta_dot, x(t), x_dot, F], [X(1, i), X(2, i), X(3, i), X(4, i), F_val]);
    x_ddot_val = subs(x_ddot, [theta(t), theta_dot, x(t), x_dot, F], [X(1, i), X(2, i), X(3, i), X(4, i), F_val]);
    X_dot = double([X(2, i); theta_ddot_val; X(4, i); x_ddot_val]);
    X = [X, X_dot * predict_step + X(:, i)];
end
J1 = 0;
Q = diag([10, 1, 10, 1]);
R = 1;
for i = 1 : 1 : predict_n + 1
    J1 = J1 + (X(:, i) - x_ref)' * Q * (X(:, i) - x_ref) + F_val' * R * F_val;
end
F_val = 10;
X = x0;
for i = 1 : 1 : predict_n
    theta_ddot_val = subs(theta_ddot, [theta(t), theta_dot, x(t), x_dot, F], [X(1, i), X(2, i), X(3, i), X(4, i), F_val]);
    x_ddot_val = subs(x_ddot, [theta(t), theta_dot, x(t), x_dot, F], [X(1, i), X(2, i), X(3, i), X(4, i), F_val]);
    X_dot = double([X(2, i); theta_ddot_val; X(4, i); x_ddot_val]);
    X = [X, X_dot * predict_step + X(:, i)];
end
J2 = 0;
Q = diag([10, 1, 10, 1]);
for i = 1 : 1 : predict_n + 1
    J2 = J2 + (X(:, i) - x_ref)' * Q * (X(:, i) - x_ref) + F_val' * R * F_val;
end
disp(J1);
disp(J2);


