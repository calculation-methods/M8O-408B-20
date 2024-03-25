import matplotlib.pyplot as plt
from numpy import *


def runge_kutta(points, start_cond_0, start_cond_1, step):
    values = [start_cond_0]
    values_der = [start_cond_1]

    for i in range(0, len(points) - 1):
        x_k = points[i]
        y_k = values[i]
        z_k = values_der[i]

        k1_1 = step * f(x_k, y_k, z_k)
        k1_2 = step * g(x_k, y_k, z_k)

        k2_1 = step * f(x_k + 0.5 * step, y_k + 0.5 * k1_1, z_k + 0.5 * k1_2)
        k2_2 = step * g(x_k + 0.5 * step, y_k + 0.5 * k1_1, z_k + 0.5 * k1_2)

        k3_1 = step * f(x_k + 0.5 * step, y_k + 0.5 * k2_1, z_k + 0.5 * k2_2)
        k3_2 = step * g(x_k + 0.5 * step, y_k + 0.5 * k2_1, z_k + 0.5 * k2_2)

        k4_1 = step * f(x_k + step, y_k + k3_1, z_k + k3_2)
        k4_2 = step * g(x_k + step, y_k + k3_1, z_k + k3_2)

        dy_k = 1 / 6 * (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1)
        dz_k = 1 / 6 * (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2)

        values.append(y_k + dy_k)
        values_der.append(z_k + dz_k)

    print(values)
    plt.plot(points, values, color='r', label='aaaaa')
    plt.show()

    return values


def get_points(begin, end, step):
    return arange(begin, end, step).tolist() + [end]


def get_values(points, func):
    return [func(points[i]) for i in range(len(points))]


def f(x, y, z):
    return -0.04 * x + 1.4 * y * z


def g(x, y, z):
    return 0.04 * x - 1.4 * y * z - 3 * y**2


def main():
    interval = [0, 10]
    step = 0.5

    points = get_points(*interval, step)
    runge_kutta(points, 1, 0.1, step)


if __name__ == '__main__':
    main()
