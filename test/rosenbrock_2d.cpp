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

  Eigen::Vector2d params(-2, -3);
  bfgs::optimizer o(rosenbrock_2d, params);
  auto x_star = o.optimize();
  std::cerr << x_star << std::endl;

  // std::cerr << "End of rosenbrock test" << std::endl;
}