import copy

import numpy as np

from functions import phi1, phi2, phi3, phi4
from show import show_inaccuracy, show_result
from task import bx, by, c, count_x, count_y, eps, hx, hy, lx, max_iterations


def liebmann_method():
    u = np.zeros((count_x + 1, count_y + 1))

    u[0, 0] = phi1(0)
    u[-1, -1] = phi2(-hy)

    for i in range(1, count_x):
        u[i, 0] = phi3(i * hx)
        u[i, -1] = phi4(i * hx)

    for j in range(1, count_y):
        u[0, j] = phi1(j * hy)
        u[-1, j] = phi2(j * hy)
        for i in range(1, count_x):
            u[i, j] = u[0, j] + (u[-1, j] - u[0, j]) / lx * i * hx

    k = 0
    while True:
        k += 1
        if k > max_iterations:
            print("Достигнуто максимальное число итераций!")
            break

        u_prev = copy.deepcopy(u)

        for j in range(1, count_y):
            for i in range(1, count_x):
                u[i, j] = (
                    -(u_prev[i + 1, j] + u_prev[i - 1, j])
                    - hx**2 * (u_prev[i, j + 1] + u_prev[i, j - 1]) / (hy**2)
                    - bx * hx * (u_prev[i + 1, j] - u_prev[i - 1, j]) / 2
                    - by * hx**2 * (u_prev[i, j + 1] - u_prev[i, j - 1]) / (2 * hy)
                ) / (c * hx**2 - 2 * (hy * hy + 1 * hx**2) / (hy**2))

        norm = np.linalg.norm(u - u_prev, np.inf)
        if norm <= eps:
            break

    print("liebmann_method: k =", k)
    return u


def sor_method(omega):
    u = np.zeros((count_x + 1, count_y + 1))

    u[0, 0] = phi1(0)
    u[-1, -1] = phi2(-hy)

    for i in range(1, count_x):
        u[i, 0] = phi3(i * hx)
        u[i, -1] = phi4(i * hx)

    for j in range(1, count_y):
        u[0, j] = phi1(j * hy)
        u[-1, j] = phi2(j * hy)
        for i in range(1, count_x):
            u[i, j] = u[0, j] + (u[-1, j] - u[0, j]) / lx * i * hx

    k = 0
    while True:
        k = k + 1
        if k > max_iterations:
            print("Достигнуто максимальное число итераций!")
            break

        u_prev = copy.deepcopy(u)

        for j in range(1, count_y):
            for i in range(1, count_x):
                u[i, j] = (
                    (
                        -(u_prev[i + 1, j] + u[i - 1, j])
                        - 1 * hx**2 * (u_prev[i, j + 1] + u[i, j - 1]) / (hy**2)
                        - bx * hx * (u_prev[i + 1, j] - u[i - 1, j]) / 2
                        - by * hx**2 * (u_prev[i, j + 1] - u[i, j - 1]) / (2 * hy)
                    )
                    / (c * hx**2 - 2 * (hy**2 + 1 * hx**2) / (hy**2))
                ) * omega + (1 - omega) * u_prev[i, j]

        norm = np.linalg.norm(u - u_prev, np.inf)
        if norm <= eps:
            break

    if omega == 1:
        print("seildel_method: k =", k)
    else:
        print("sor_method: k =", k)

    return u


def get_axis_np(count, mul):
    axis = np.zeros(count)
    for i in range(count):
        axis[i] = mul * i
    return axis


def main():
    u1 = liebmann_method()
    u2 = sor_method(1)
    u3 = sor_method(1.5)

    y_axis = get_axis_np(count_y + 1, hy)
    x_axis = get_axis_np(count_x + 1, hx)

    show_result(y_axis, x_axis, u1, u2, u3)
    show_inaccuracy(y_axis, x_axis, u1)


if __name__ == "__main__":
    main()
