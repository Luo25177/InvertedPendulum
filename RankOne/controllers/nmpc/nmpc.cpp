#include "inverted_pendulum.hpp"
#include <cstdio>

int main(int argc, char **argv) {
  InvertedPendulum *robot = new InvertedPendulum();

  FILE *fp;
  fopen_s(&fp, "../../data/Data.xlsx", "w");

  while (robot->step(time_step) != -1) {
    robot->update();
    robot->show(fp);
    robot->run_nmpc();
  };

  fclose(fp);
  delete robot;
  return 0;
}
