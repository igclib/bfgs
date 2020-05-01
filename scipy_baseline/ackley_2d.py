from scipy.optimize import minimize
import numpy as np


def ackley_2d(x):
    return np.sum(100.0 * (x[1:] - x[:-1] ** 2.0) ** 2.0 + (1 - x[:-1]) ** 2.0, axis=0)


if __name__ == "__main__":
    init_vector = [0.6, -0.4]
    res = minimize(ackley_2d, init_vector, tol=1e-5)
    print(res)
