#include <bfgs>
#include <iostream>

double rosenbrock_2d(double *x);
void plot();

void plot() {
  for (double i = -3; i < 3; i += 0.1) {
    for (double j = -3; j < 3; j += 0.1) {
      double params[2] = {i, j};
      double res = rosenbrock_2d(params);
      // std::cout << "[" << i << "," << j << "] " << res << std::endl;
      std::cout << i << " " << j << " " << res << std::endl;
    }
    std::cout << std::endl;
  }
}

double rosenbrock_2d(double *x) {
  return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) +
         (1 - x[0]) * (1 - x[0]);
}

int main() { plot(); }