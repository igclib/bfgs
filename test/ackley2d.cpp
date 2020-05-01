#include <bfgs.hpp>
#include <eigen3/Eigen/Core>
#include <fstream>
#include <iostream>
#include <vector>

double ackley2d(const Eigen::VectorXd &x);
void plot_f();
void plot_dist();
void plot_path();

double ackley2d(const Eigen::VectorXd &x) {
  return -200 * std::exp(-0.02 * std::sqrt((x(0) * x(0)) + (x(1) * x(1)))) +
         5 * std::exp(std::cos(3 * x(0)) + std::sin(3 * x(1)));
}

void plot_f() {
  std::ofstream f;
  f.open("plot/ackley2d");
  for (double i = -4; i <= 4; i += 0.1) {
    for (double j = -4; j <= 4; j += 0.1) {
      Eigen::Vector2d params(i, j);
      double res = ackley2d(params);
      f << i << " " << j << " " << res << std::endl;
    }
    f << std::endl;
  }
}

void plot_dist() {
  std::ofstream f;
  f.open("plot/ackleydist");
  for (double i = -5; i <= 5; i += 0.1) {
    for (double j = -5; j <= 5; j += 0.1) {
      auto params = Eigen::Vector2d(i, j);
      bfgs::optimizer o(ackley2d, params);
      auto result = o.optimize();
      f << result.x(0) << ' ' << result.x(1) << std::endl;
    }
  }
}

void plot_path() {
  std::ofstream f;
  f.open("plot/ackleypath");
  auto params = Eigen::Vector2d(2, 2);
  bfgs::optimizer o(ackley2d, params, true);
  auto result = o.optimize();
  for (auto x : result.path) {
    f << x(0) << " " << x(1) << " " << ackley2d(x) << std::endl;
  }
}

int main() {

  plot_f();
  plot_dist();
  plot_path();
}