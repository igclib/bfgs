from scipy.optimize import minimize
import numpy as np


def rosenbrok_2d(x):
    return np.sum(100.0 * (x[1:] - x[:-1] ** 2.0) ** 2.0 + (1 - x[:-1]) ** 2.0, axis=0)


if __name__ == "__main__":
    init_vector = [-3, -0.89999999999999836]
    res = minimize(rosenbrok_2d, init_vector, tol=1e-5)
    print(res)
