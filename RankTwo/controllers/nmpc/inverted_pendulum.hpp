#pragma once

#include "nmpc_solver.hh"
#include <webots/Motor.hpp>
#include <webots/PositionSensor.hpp>
#include <webots/Robot.hpp>

#define time_step 5
#define delta_time 0.005

class InvertedPendulum : public webots::Robot {
  webots::Motor*          motor;
  webots::PositionSensor *x_sensor, *theta_sensor, *phi_sensor;

  double x;
  double x_dot;
  double theta;
  double theta_dot;
  double phi;
  double phi_dot;
  double f;

  NMPCSolver nmpc_solver;

public:
  InvertedPendulum();

  ~InvertedPendulum() { }

  void update();
  void run_lqr();
  void run_nmpc();
};
