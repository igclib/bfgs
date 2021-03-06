#include <bfgs.hpp>
#include <eigen3/Eigen/Core>
#include <fstream>
#include <iostream>
#include <vector>

double rosenbrock_2d(const Eigen::VectorXd &x);
void plot_f();
void plot_dist();
void plot_path();

void plot_f() {
  std::ofstream f;
  f.open("plot/rosenbrock2d");
  for (double i = -5; i <= 5; i += 0.1) {
    for (double j = -5; j <= 5; j += 0.1) {
      Eigen::Vector2d params(i, j);
      double res = rosenbrock_2d(params);
      f << i << " " << j << " " << res << std::endl;
    }
    f << std::endl;
  }
}

void plot_dist() {
  std::ofstream f;
  f.open("plot/rosdist");
  for (double i = -5; i <= 5; i += 0.1) {
    for (double j = -5; j <= 5; j += 0.1) {
      auto params = Eigen::Vector2d(i, j);
      bfgs::optimizer o(rosenbrock_2d, params);
      auto result = o.optimize();
      f << result.x(0) << ' ' << result.x(1) << std::endl;
    }
  }
}

void plot_path() {
  std::ofstream f;
  f.open("plot/rospath");
  auto params = Eigen::Vector2d(-2, -2);
  bfgs::optimizer o(rosenbrock_2d, params, true);
  auto result = o.optimize();
  for (auto x : result.path) {
    f << x(0) << " " << x(1) << " " << rosenbrock_2d(x) << std::endl;
  }
}

double rosenbrock_2d(const Eigen::VectorXd &x) {
  return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) +
         (1 - x[0]) * (1 - x[0]);
}

int main() {

  plot_f();
  plot_dist();
  plot_path();
}