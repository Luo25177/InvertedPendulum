#include "inverted_pendulum.hpp"

using namespace webots;

int main(int argc, char** argv) {
  InvertedPendulum* robot = new InvertedPendulum();

  while (robot->step(time_step) != -1) {
    robot->update();
    robot->run_nmpc();
  };

  delete robot;
  return 0;
}
