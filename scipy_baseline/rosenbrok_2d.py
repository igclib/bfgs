from scipy.optimize import minimize


def rosenbrok_2d(x):
    return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0])


if __name__ == "__main__":
    init_vector = [0, 0]
    res = minimize(rosenbrok_2d, init_vector, tol=1e-5)
    print(res)
