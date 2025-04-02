#include "inverted_pendulum.hpp"

#include <Eigen>

#define nmpc_h 30
#define nmpc_step 0.03

InvertedPendulum::InvertedPendulum() {
  motor = new webots::Motor("BasePosMotor");
  x_sensor = new webots::PositionSensor("BasePosSensor");
  theta_sensor = new webots::PositionSensor("Rob1Sensor");
  phi_sensor = new webots::PositionSensor("Rob2Sensor");

  x_sensor->enable(time_step);
  theta_sensor->enable(time_step);
  phi_sensor->enable(time_step);

  DM Q = DM::zeros(6, 6);
  DM R = DM::zeros(1, 1);

  Q(0, 0) = 1000;
  Q(1, 1) = 100;
  Q(2, 2) = 1000;
  Q(3, 3) = 100;
  Q(4, 4) = 1000;
  Q(5, 5) = 100;

  R(0, 0) = 0.00005;

  SX X = SX::sym("x", 6);
  SX U = SX::sym("u", 1);

  SX theta = X(0);
  SX theta_dot = X(1);
  SX x = X(2);
  SX x_dot = X(3);
  SX phi = X(4);
  SX phi_dot = X(5);
  SX F = U(0);

  // clang-format off
  SX numerator_theta = -(1./*441 * (
      706.3 * sin(theta) +
      429.6 * pow(cos(phi), 2) * sin(theta) -
      3.459 * pow(theta_dot, 2) * sin(2.0 * phi) -
      2.679 * pow(theta_dot, 2) * sin(2.0 * theta) -
      37.04 * F * cos(theta) -
      429.6 * cos(phi) * cos(theta) * sin(phi) +
      18.45 * F * pow(cos(phi), 2) * cos(theta) -
      15.45 * pow(phi_dot, 2) * cos(phi) * sin(theta) +
      11.73 * pow(phi_dot, 2) * cos(theta) * sin(phi) +
      24.44 * pow(theta_dot, 2) * cos(phi) * pow(cos(theta), 2) * sin(phi) -
      26.29 * pow(theta_dot, 2) * pow(cos(phi), 2) * cos(theta) * sin(theta) +
      18.45 * F * cos(phi) * sin(phi) * sin(theta)));
  SX denominator_theta = (
      7.723 * pow(cos(theta), 2) -
      27.01 * pow(cos(phi), 2) +
      37.88 * pow(cos(phi), 2) * pow(cos(theta), 2) +
      35.22 * cos(phi) * cos(theta) * sin(phi) * sin(theta) -
      32.41);
  SX theta_ddot = numerator_theta / denominator_theta;

  SX numerator_x = -(0.32 * (
        26.13 * F +
        26.94 * sin(2.0 * phi) -
        164.2 * sin(2.0 * theta) +
        33.23 * F * pow(cos(phi), 2) +
        10.99 * F * pow(cos(theta), 2) +
        5.227 * pow(phi_dot, 2) * sin(phi) +
        10.45 * pow(theta_dot, 2) * sin(theta) -
        6.698 * pow(phi_dot, 2) * pow(cos(theta), 2) * sin(phi) +
        4.422 * pow(theta_dot, 2) * pow(cos(phi), 2) * sin(theta) -
        81.46 * pow(cos(phi), 2) * cos(theta) * sin(theta) -
        49.85 * F * pow(cos(phi), 2) * pow(cos(theta), 2) +
        5.573 * pow(phi_dot, 2) * cos(phi) * cos(theta) * sin(theta) -
        2.198 * pow(theta_dot, 2) * cos(phi) * cos(theta) * sin(phi) -
        4.985 * F * cos(phi) * cos(theta) * sin(phi) * sin(theta)));

  SX denominator_x = (
      7.723 * pow(cos(theta), 2) -
      27.01 * pow(cos(phi), 2) +
      37.88 * pow(cos(phi), 2) * pow(cos(theta), 2) +
      35.22 * cos(phi) * cos(theta) * sin(phi) * sin(theta) -
      32.41);

  SX x_ddot = numerator_x / denominator_x;

  SX numerator_phi = 1.441 * (
      231.8 * sin(phi) -
      9.372 * pow(phi_dot, 2) * sin(2.0 * phi) +
      8.762 * pow(phi_dot, 2) * sin(2.0 * theta) -
      859.1 * pow(cos(theta), 2) * sin(phi) +
      49.24 * F * cos(phi) +
      508.7 * cos(phi) * cos(theta) * sin(theta) -
      46.12 * F * cos(phi) * pow(cos(theta), 2) -
      18.46 * pow(theta_dot, 2) * cos(phi) * sin(theta) +
      25.59 * pow(theta_dot, 2) * cos(theta) * sin(phi) +
      26.29 * pow(phi_dot, 2) * cos(phi) * pow(cos(theta), 2) * sin(phi) -
      24.44 * pow(phi_dot, 2) * pow(cos(phi), 2) * cos(theta) * sin(theta) -
      55.34 * F * cos(theta) * sin(phi) * sin(theta)
  );

  SX denominator_phi = (
      7.723 * pow(cos(theta), 2) -
      27.01 * pow(cos(phi), 2) +
      37.88 * pow(cos(phi), 2) * pow(cos(theta), 2) +
      35.22 * cos(phi) * cos(theta) * sin(phi) * sin(theta) -
      32.41
  );

  SX phi_ddot =  numerator_phi / denominator_phi;

  SX X_dot = vertcat(
    theta_dot,
    theta_ddot,
    x_dot,
    x_ddot,
    phi_dot,
    phi_ddot*/
  //);
  // clang-format on
  
  DM A = DM::zeros(6, 6); 
  A(0, 1) = 1;
  A(2, 3) = 1;
  A(4, 5) = 1;
  A(1, 0) = 118.494123950081;
  A(3, 0) = -9.49279645740227;
  A(5, 0) = -53.0653165379936;
  A(1, 4) = -44.8107117431946;
  A(3, 4) = 1.24801747245928;
  A(5, 4) = 65.4410849574907;

  DM B = DM::zeros(6, 1);

  B(1, 0) = -1.93932423436066;
  B(3, 0) = 0.475064505862854;
  B(5, 0) = -0.325726332371490;

  SX X_dot = SX::mtimes(A, X) + SX::mtimes(B, U);

  Function dynamics = Function("pendulum", vector<SX> { X, U }, vector<SX> { X_dot }, vector<string> { "x", "u" }, vector<string> { "xdot" });

  nmpc_solver.setup(nmpc_h, nmpc_step, 6, 1, dynamics, Q, R, -20, 20);
}

void InvertedPendulum::update() {
  double last_x = x;
  double last_theta = theta;
  double last_phi = phi;

  x = x_sensor->getValue();
  theta = theta_sensor->getValue();
  phi = phi_sensor->getValue() + theta;

  x_dot = (x - last_x) / delta_time;
  theta_dot = (theta - last_theta) / delta_time;
  phi_dot = (phi - last_phi) / delta_time;
}

void InvertedPendulum::run_lqr() {
  Eigen::Matrix<double, 6, 1> X;
  Eigen::Matrix<double, 6, 1> X_ref;
  Eigen::Matrix<double, 1, 1> U;
  Eigen::Matrix<double, 1, 6> K;

  X << theta, theta_dot, x, x_dot, phi, phi_dot;
  X_ref << 0, 0, 0.5, 0, 0, 0;
  K << -589.445674284356, -32.5650180786876, 10.0000000000010, 17.3792120847860, 694.622444754623, 80.6098143116669;

  U = K * (X_ref - X);

  motor->setForce(U(0, 0));
}

void InvertedPendulum::run_nmpc() {
  DM x0 = DM({ theta, theta_dot, x, x_dot, phi, phi_dot });
  DM x_ref = DM({ 0, 0, 0.5, 0, 0, 0 });
  DM u = nmpc_solver.solve(x0, x_ref);
  f = static_cast<double>(u(0, 0));
  motor->setForce(f);
}
