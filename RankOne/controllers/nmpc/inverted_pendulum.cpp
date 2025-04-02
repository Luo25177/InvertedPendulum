#include "inverted_pendulum.hpp"

#define nmpc_h 20
#define nmpc_step 0.025

InvertedPendulum::InvertedPendulum() {
  motor = new webots::Motor("BasePosMotor");
  x_sensor = new webots::PositionSensor("BasePosSensor");
  theta_sensor = new webots::PositionSensor("Rob1Sensor");

  x_sensor->enable(time_step);
  theta_sensor->enable(time_step);
  DM Q = DM::zeros(4, 4);
  DM R = DM::zeros(1, 1);

  Q(0, 0) = 100;
  Q(1, 1) = 10;
  Q(2, 2) = 100;
  Q(3, 3) = 10;

  R(0, 0) = 0.02;

  SX X = SX::sym("x", 4);
  SX U = SX::sym("u", 1);

  SX theta = X(0);
  SX theta_dot = X(1);
  SX x = X(2);
  SX x_dot = X(3);
  SX F = U(0);

  SX numerator_theta = 0.95 * cos(theta) * sin(theta) * theta_dot * theta_dot + 28.1875 * sin(theta) - F * cos(theta);
  SX denominator_theta = cos(theta) * cos(theta) - 0.4004;
  SX theta_ddot = 1.0524 * (numerator_theta / denominator_theta);

  SX numerator_x = 1413.0 * cos(theta) * sin(theta) -48.79 * F + 144.1 * F * (pow(sin(theta), 2) - pow(cos(theta), 2)) + 19.062 * pow(theta_dot, 2) * sin(theta);
  SX denominator_x = 13.136 * pow(cos(theta), 2) - 8.771 * pow(sin(theta), 2);
  SX x_ddot = -0.032 * numerator_x / denominator_x;

  SX X_dot = vertcat(
    theta_dot,
    theta_ddot,
    x_dot,
    x_ddot
  );

  Function dynamics = Function("pendulum", std::vector<SX>{X, U}, std::vector<SX>{X_dot}, std::vector<std::string>{"x", "u"}, std::vector<std::string>{"xdot"});

  nmpc_solver.setup(nmpc_h, nmpc_step, 4, 1, dynamics, Q, R, -20, 20);
}


void InvertedPendulum::update() {
  double last_x = x;
  double last_theta = theta;

  x = x_sensor->getValue();
  theta = theta_sensor->getValue();

  x_dot = (x - last_x) / delta_time;
  theta_dot = (theta - last_theta) / delta_time;
}

void InvertedPendulum::run_lqr() {
  Eigen::Matrix<double, 4, 1> X;
  Eigen::Matrix<double, 4, 1> X_ref;
  Eigen::Matrix<double, 1, 1> U;
  Eigen::Matrix<double, 1, 4> K;

  X << theta, theta_dot, x, x_dot;
  X_ref << 0, 0, 0.3, 0;
  K << -84.0876402902565, -12.5574101757239, -9.99999999999990, -10.7247339827630;

  U = K * (X_ref - X);
  motor->setForce(U(0, 0));
}

void InvertedPendulum::run_nmpc() {
  DM x0 = DM({theta, theta_dot, x, x_dot});
  DM x_ref = DM({ 0, 0, 0.5, 0 });
  DM u = nmpc_solver.solve(x0, x_ref);
  f = static_cast<double>(u(0, 0));
  motor->setForce(f);
}

void InvertedPendulum::show(FILE* fp) {
  fprintf(fp, "%f\t%f\t%f\t%f\n", theta, theta_dot, x, x_dot);
}

