import numpy as np

from functions import f, phi1, phi2, phi3, phi4, psi
from show import show_inaccuracy, show_result
from task import count_x, count_y, hx, hy, t_count, tau


def alternative_directions_scheme():
    u = np.zeros((count_x + 1, count_y + 1, t_count + 1))
    u_1 = np.zeros((count_x + 1, count_y + 1))
    u_2 = np.zeros((count_x + 1, count_y + 1))

    ai = np.zeros(count_x + 1)
    bi = np.zeros(count_x + 1)
    ci = np.zeros(count_x + 1)
    di = np.zeros(count_x + 1)

    for i in range(count_x + 1):
        for j in range(count_y + 1):
            u[i, j, 0] = psi(i * hx, j * hy)

    for k in range(1, t_count + 1):
        u_prev = u[:, :, k - 1]
        t_step = tau * (k - 0.5)

        for j in range(count_y):
            bi[0] = hx
            bi[-1] = hx - 1
            ci[0] = 0
            ai[-1] = 1
            di[0] = phi1(j * hy, t_step) * hx
            di[-1] = phi2(j * hy, t_step) * hx
            for i in range(1, count_x):
                ai[i] = 1
                bi[i] = -2 * (hx**2) / tau - 2
                ci[i] = 1
                di[i] = (
                    -2 * (hx**2) * u_prev[i, j] / tau
                    - (hx**2)
                    * (u_prev[i, j + 1] - 2 * u_prev[i, j] + u_prev[i, j - 1])
                    / (hy**2)
                    - (hx**2) * f(i * hx, j * hy, t_step)
                )

            ta = thomas_algorithm(ai, bi, ci, di)
            for i in range(count_x + 1):
                u_1[i, j] = ta[i]
                u_1[i, 0] = phi3(i * hx, t_step)
                u_1[i, -1] = (phi4(i * hx, t_step) - u_1[i, -2] / hy) / (1 - 1 / hy)

        for j in range(count_y + 1):
            u_1[0, j] = phi1(j * hy, t_step)
            u_1[-1, j] = (phi2(j * hy, t_step) - u_1[-2, j] / hx) / (1 - 1 / hx)

        for i in range(count_x):
            bi[0] = hy
            bi[-1] = hy - 1
            ci[0] = 0
            ai[-1] = 1
            di[0] = phi3(i * hx, k * tau) * hy
            di[-1] = phi4(i * hx, k * tau) * hy

            for j in range(1, count_y):
                ai[j] = 1
                bi[j] = -2 * (hy**2) / tau - 2
                ci[j] = 1
                di[j] = (
                    -2 * (hy**2) * u_1[i, j] / tau
                    - (hy**2)
                    * (u_1[i + 1, j] - 2 * u_1[i, j] + u_1[i - 1, j])
                    / (hx**2)
                    - (hy**2) * f(i * hx, j * hy, k * tau)
                )

            ta = thomas_algorithm(ai, bi, ci, di)
            for j in range(count_y + 1):
                u_2[i, j] = ta[j]
                u_2[0, j] = phi1(j * hy, k * tau)
                u_2[-1, j] = (phi2(j * hy, k * tau) - u_2[-2, j] / hx) / (1 - 1 / hx)

        for i in range(count_x + 1):
            u_2[i, 0] = phi3(i * hx, k * tau)
            u_2[i, -1] = (phi4(i * hx, k * tau) - u_2[i, -2] / hy) / (1 - 1 / hy)
            for j in range(count_y + 1):
                u[i, j, k] = u_2[i, j]

    return u


# схема дробных шагов
def fractional_steps_scheme():
    u = np.zeros((count_x + 1, count_y + 1, t_count + 1))
    u_1 = np.zeros((count_x + 1, count_y + 1))
    u_2 = np.zeros((count_x + 1, count_y + 1))

    ai = np.zeros(count_x + 1)
    bi = np.zeros(count_x + 1)
    ci = np.zeros(count_x + 1)
    di = np.zeros(count_x + 1)

    for i in range(count_x + 1):
        for j in range(count_y + 1):
            u[i, j, 0] = psi(i * hx, j * hy)

    for k in range(1, t_count + 1):
        u_prev = u[:, :, k - 1]
        t_step = tau * (k - 1)

        for j in range(count_y):
            bi[0] = hx
            bi[-1] = hx - 1
            ci[0] = 0
            ai[-1] = 1
            di[0] = phi1(j * hy, t_step) * hx
            di[-1] = phi2(j * hy, t_step) * hx
            for i in range(1, count_x):
                ai[i] = 1
                bi[i] = -(hx**2) / tau - 2
                ci[i] = 1
                di[i] = (
                    -(hx**2) * u_prev[i, j] / tau
                    - (hx**2) * f(i * hx, j * hy, t_step) / 2
                )

            ta = thomas_algorithm(ai, bi, ci, di)
            for i in range(count_x + 1):
                u_1[i, j] = ta[i]
                u_1[i, 0] = phi3(i * hx, t_step)
                u_1[i, -1] = (phi4(i * hx, t_step) - u_1[i, -2] / hy) / (1 - 1 / hy)

        for j in range(count_y + 1):
            u_1[0, j] = phi1(j * hy, t_step)
            u_1[-1, j] = (phi2(j * hy, t_step) - u_1[-2, j] / hx) / (1 - 1 / hx)

        for i in range(count_x):
            bi[0] = hy
            bi[-1] = hy - 1
            ci[0] = 0
            ai[-1] = 1
            di[0] = phi3(i * hx, k * tau) * hy
            di[-1] = phi4(i * hx, k * tau) * hy

            for j in range(1, count_y):
                ai[j] = 1
                bi[j] = -(hy**2) / tau - 2
                ci[j] = 1
                di[j] = (
                    -(hy**2) * u_1[i, j] / tau
                    - (hy**2) * f(i * hx, j * hy, k * tau) / 2
                )
            ta = thomas_algorithm(ai, bi, ci, di)
            for j in range(count_y + 1):
                u_2[i, j] = ta[j]
                u_2[0, j] = phi1(j * hy, k * tau)
                u_2[-1, j] = (phi2(j * hy, k * tau) - u_2[-2, j] / hx) / (1 - 1 / hx)

        for i in range(count_x + 1):
            u_2[i, 0] = phi3(i * hx, k * tau)
            u_2[i, -1] = (phi4(i * hx, k * tau) - u_2[i, -2] / hy) / (1 - 1 / hy)
            for j in range(count_y + 1):
                u[i, j, k] = u_2[i, j]

    return u


def thomas_algorithm(a, b, c, d):
    size = len(a)
    P = np.zeros(size)
    Q = np.zeros(size)
    P[0] = -c[0] / b[0]
    Q[0] = d[0] / b[0]

    for i in range(1, size):
        s = b[i] + a[i] * P[i - 1]
        P[i] = -c[i] / s
        Q[i] = (d[i] - a[i] * Q[i - 1]) / s

    result = np.zeros(size)
    result[-1] = Q[-1]

    for i in range(size - 2, -1, -1):
        result[i] = P[i] * result[i + 1] + Q[i]

    return result


#
def main():
    u1 = alternative_directions_scheme()
    u2 = fractional_steps_scheme()

    show_result(u1, u2)
    show_inaccuracy(u1, u2)


if __name__ == "__main__":
    main()
