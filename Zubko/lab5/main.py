import numpy as np

from functions import phi0, phi1, psi
from show import show_inaccuracy, show_result
from task import a, c, count_t, count_x, h, sigma, tau


def explicit_scheme(bound_condition):
    u = np.zeros((count_t, count_x))

    for i in range(1, count_x):
        u[0][i] = psi(i * h)

    for k in range(1, count_t):
        for i in range(1, count_x - 1):
            u[k][i] = (
                sigma * u[k - 1][i + 1]
                + (1 - 2 * sigma) * u[k - 1][i]
                + sigma * u[k - 1][i - 1]
                + c * tau * u[k - 1][i]
            )

        if bound_condition == 1:
            u[k][0] = u[k][1] - h * phi0(k * tau)
            u[k][-1] = phi1(k * tau)
        elif bound_condition == 2:
            u[k][0] = (phi0(k * tau) + u[k][2] / (2 * h) - 2 * u[k][1] / h) * 2 * h / -3
            u[k][-1] = phi1(k * tau)
        elif bound_condition == 3:
            u[k][0] = (
                u[k][1] - h * phi0(k * tau) + (h**2 / (2 * tau) * u[k - 1][0])
            ) / (1 + h**2 / (2 * tau))
            u[k][-1] = phi1(k * tau)
        else:
            print("Условие не найдено")
    return u


def implicit_scheme(bound_condition):
    u = np.zeros((count_t, count_x))

    ai = np.zeros(count_x)
    bi = np.zeros(count_x)
    ci = np.zeros(count_x)
    di = np.zeros(count_x)

    for i in range(1, count_x):
        u[0][i] = psi(i * h)

    for k in range(1, count_t):
        for i in range(1, count_x - 1):
            ai[i] = sigma
            bi[i] = -2 * sigma + c * tau - 1
            ci[i] = sigma
            di[i] = -u[k - 1][i]

        if bound_condition == 1:
            bi[0] = -(1 + 2 * sigma - c * tau)
            ci[0] = 2 * sigma
            di[0] = -(u[k - 1][0] - 2 * a * tau * phi0(k * tau) / h)
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


def explicit_implicit_scheme(omega, bound_condition):
    u = np.zeros((count_t, count_x))

    imp = implicit_scheme(bound_condition)
    exp = explicit_scheme(bound_condition)

    for k in range(count_t):
        for i in range(count_x):
            u[k][i] = imp[k][i] * omega + exp[k][i] * (1 - omega)

    return u


def get_axis_np(count, mul):
    axis = np.zeros(count)
    for i in range(count):
        axis[i] = mul * i
    return axis


def solve():
    res1 = explicit_scheme(1)
    res2 = implicit_scheme(1)
    res3 = explicit_implicit_scheme(0.5, 1)

    t_axis = get_axis_np(count_t, tau)
    x_axis = get_axis_np(count_x, h)

    show_result(t_axis, x_axis, res1, res2, res3)
    show_inaccuracy(t_axis, x_axis, res1)


if __name__ == "__main__":
    solve()
