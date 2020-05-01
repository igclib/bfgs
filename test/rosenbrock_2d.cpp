#include <bfgs>
#include <eigen3/Eigen/Core>
#include <iostream>
#include <vector>

double rosenbrock_2d(const Eigen::VectorXd &x);
void plot();

void plot() {
  for (double i = -3; i < 3; i += 0.1) {
    for (double j = -3; j < 3; j += 0.1) {
      Eigen::Vector2d params(i, j);
      double res = rosenbrock_2d(params);
      std::cout << i << " " << j << " " << res << std::endl;
    }
    std::cout << std::endl;
  }
}

double rosenbrock_2d(const Eigen::VectorXd &x) {
  return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) +
         (1 - x[0]) * (1 - x[0]);
}

int main() {

  // plot();
  /*
    auto params = Eigen::Vector2d(-3, -0.89999999999999836);
    bfgs::optimizer o(rosenbrock_2d, params);
    auto result = o.optimize();
    std::cerr << result << std::endl;
  */

  const double MAXTOL = 1e-4;
  Eigen::Vector2d expected(1.0, 1.0);
  for (double i = -3; i < 3; i += 0.1) {
    for (double j = -3; j < 3; j += 0.1) {
      auto params = Eigen::Vector2d(i, j);
      bfgs::optimizer o(rosenbrock_2d, params);
      auto result = o.optimize();
      double error = (expected - result.x).cwiseAbs().maxCoeff();
      if (!result.success) {
        std::cerr << "No convergence for " << i << ", " << j << std::endl;
      }
      if (error > MAXTOL) {
        std::cerr << error << std::endl;
        std::cerr << "Too far for " << i << ", " << j << std::endl;
      }
    }
  }
}