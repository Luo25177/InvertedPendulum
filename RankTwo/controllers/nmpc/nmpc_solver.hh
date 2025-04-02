#pragma once

#include <casadi/casadi.hpp>

#include <map>
#include <string>
#include <vector>

using namespace std;
using namespace casadi;

class NMPCSolver {
  int    state_n;
  int    control_n;
  int    predict_n;
  double predict_step;

  Function        solver;
  map<string, DM> args;
  vector<double>  res;

public:
  NMPCSolver();
  ~NMPCSolver();

  void setup(const int &_predict_n, const double &_predict_step, const int &_state_n, const int &_control_n, Function _dynamics, const DM &_Q, const DM &_R, const DM &_u_min, const DM &_u_max);

  DM solve(const DM &x0, const DM &x_ref);
};
