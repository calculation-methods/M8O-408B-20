import numpy as np

from functions import phi0, phi1, psi1, psi2
from show import show_inaccuracy, show_result
from task import a, b, c, count_t, count_x, e, h, sigma, tau


def explicit_scheme(bound_condition, initial_condition):
    u = np.zeros((count_t, count_x))

    for i in range(0, count_x):
        u[0][i] = psi1(i * h)

        if initial_condition == 1:
            u[1][i] = psi1(i * h) + psi2(i * h) * tau + 0
        elif initial_condition == 2:
            u[1][i] = (
                psi1(i * h)
                + tau * psi2(i * h)
                + 0.5 * (tau**2) * (-e * psi2(i * h) + c * psi1(i * h))
            )
        else:
            print("Approximation type not found")

    for k in range(2, count_t):
        for i in range(1, count_x - 1):
            u[k][i] = (
                4 * u[k - 1][i]
                + (e * tau - 2) * u[k - 2][i]
                + 2
                * a
                * (tau**2)
                * (u[k - 1][i - 1] - 2 * u[k - 1, i] + u[k - 1][i + 1])
                / (h**2)
                + b * (tau**2) * (u[k - 1][i + 1] - u[k - 1][i - 1]) / h
                + 2 * (tau**2) * (c * u[k - 1][i])
            ) / (2 + e * tau)
        if bound_condition == 1:
            u[k][0] = 0
            u[k][-1] = 0
        elif bound_condition == 2:
            u[k][0] = (phi0(k * tau) + u[k][2] / (2 * h) - 2 * u[k][1] / h) * 2 * h / -3
            u[k][-1] = phi1(k * tau)
        elif bound_condition == 3:
            u[k][0] = 2 * h * phi1(k * tau)
            u[k][-1] = 2 * h * phi1(k * tau)
        else:
            print("Условие не найдено")

    return u


def implicit_scheme(bound_condition, initial_condition):
    u = np.zeros((count_t, count_x))

    ai = np.zeros(count_x)
    bi = np.zeros(count_x)
    ci = np.zeros(count_x)
    di = np.zeros(count_x)

    for i in range(0, count_x):
        u[0][i] = psi1(i * h)

        if initial_condition == 1:
            u[1][i] = psi1(i * h) + psi2(i * h) * tau
        elif initial_condition == 2:
            u[1][i] = (
                psi1(i * h)
                + tau * psi2(i * h)
                + 0.5 * (tau**2) * (-e * psi2(i * h) + c * psi1(i * h))
            )
        else:
            print("Условие не найдено")

    for k in range(2, count_t):
        for i in range(1, count_x - 1):
            ai[i] = 2 * a - h * b
            bi[i] = 2 * (h**2) * (-e / (2 * tau) - 1 / (tau**2) + c) - 4 * a
            ci[i] = h * b + 2 * a
            di[i] = (
                -4 * (h**2) * u[k - 1][i] / (tau**2)
                + (2 * (h**2) / (tau**2) - e * (h**2) / tau) * u[k - 2][i]
            )

        if bound_condition == 1:
            bi[0] = h
            ci[0] = 0
            di[0] = h * phi0(k * tau)
            ai[-1] = 2 * sigma
            bi[-1] = -(1 + 2 * sigma - c * tau)
            di[-1] = -phi1(k * tau)
        elif bound_condition == 2:
            bi[0] = -(1 + 2 * sigma - c * tau)
            ci[0] = 2 * sigma
            di[0] = -(u[k - 1][0] - 2 * a * tau * phi0(k * tau) / h)
            ai[-1] = 2 * sigma
            bi[-1] = -(1 + 2 * sigma - c * tau)
            di[-1] = -phi1(k * tau)
        elif bound_condition == 3:
            bi[0] = -(1 + 2 * sigma - c * tau)
            ci[0] = 2 * sigma
            di[0] = -(
                (1 - sigma) * u[k - 1][1] + sigma / 2 * u[k - 1][0]
            ) - sigma * phi0(k * tau)
            ai[-1] = 2 * sigma
            bi[-1] = -(1 + 2 * sigma - c * tau)
            di[-1] = -phi1(k * tau)
        else:
            print("Условие не найдено")

        u[k] = thomas_algorithm(ai, bi, ci, di)

    return u


def thomas_algorithm(a, b, c, d):
    size = len(a)
    p = np.zeros(size)
    q = np.zeros(size)
    p[0] = -c[0] / b[0]
    q[0] = d[0] / b[0]

    for i in range(1, size):
        s = b[i] + a[i] * p[i - 1]
        p[i] = -c[i] / s
        q[i] = (d[i] - a[i] * q[i - 1]) / s

    result = np.zeros(size)
    result[-1] = q[-1]

    for i in range(size - 2, -1, -1):
        result[i] = p[i] * result[i + 1] + q[i]

    return result


def get_axis_np(count, mul):
    axis = np.zeros(count)
    for i in range(count):
        axis[i] = mul * i
    return axis


def main():
    res1 = explicit_scheme(1, 1)
    res2 = implicit_scheme(1, 1)

    t_axis = get_axis_np(count_t, tau)
    x_axis = get_axis_np(count_x, h)

    show_result(t_axis, x_axis, res1, res2)
    show_inaccuracy(t_axis, x_axis, res1)


if __name__ == "__main__":
    main()
