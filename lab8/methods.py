from parameters import *


def tridiagonal_matrix_algorithm(aa, bb, cc, dd):
    dl = len(dd)
    x_ans = np.zeros(dl)
    p = [((-1)*cc[0]) / bb[0]]
    q = [dd[0] / bb[0]]

    for i in range(1, dl):
        p.append(((-1)*cc[i]) / (bb[i] + aa[i] * p[i - 1]))
        q.append((dd[i] - aa[i] * q[i-1]) / (bb[i] + aa[i] * p[i-1]))
    x_ans[-1] = q[-1]

    for i in reversed(range(dl - 1)):
        x_ans[i] = p[i] * x_ans[i+1] + q[i]

    return x_ans


def mdn():
    u = np.zeros((lx, ly, lt))

    for i in range(lx):
        for j in range(ly):
            u[i][j][0] = psi()

    for k in range(1, lt):
        u1 = np.zeros((lx, ly))
        t_half = t[k] - tau / 2
        u_k_1 = u[:, :, k - 1]

        for j in range(ly - 1):
            aa, bb, cc, dd = np.zeros(lx), np.zeros(lx), np.zeros(lx), np.zeros(lx)

            bb[0] = hx * alpha_2 - alpha_1
            cc[0] = alpha_1
            dd[0] = phi_1() * hx

            aa[-1] = - beta_1
            bb[-1] = hx * beta_2 + beta_1
            dd[-1] = phi_2(y[j], t_half) * hx

            for i in range(1, lx - 1):
                aa[i] = a
                bb[i] = -2 * (hx ** 2) / tau - 2 * a
                cc[i] = a

                dd[i] = - b * (hx ** 2) * (u_k_1[i][j+1] - 2 * u_k_1[i][j] + u_k_1[i][j-1]) / (hy ** 2) \
                        - 2 * (hx ** 2) * u_k_1[i][j] / tau \
                        - (hx ** 2) * f(x[i], y[j], t_half)

            xx = tridiagonal_matrix_algorithm(aa, bb, cc, dd)

            for i in range(lx):
                u1[i][j] = xx[i]
                u1[i][0] = (phi_3() - alpha_3 * u1[i][1] / hy) / (beta_3 - alpha_3 / hy)
                u1[i][-1] = (phi_4(x[i], t_half) + alpha_4 * u1[i][-2] / hy) / (beta_4 + alpha_4 / hy)

        for j in range(ly):
            u1[0][j] = (phi_1() - alpha_1 * u1[1][j] / hx) / (alpha_2 - alpha_1 / hx)
            u1[-1][j] = (phi_2(y[j], t_half) + beta_1 * u1[-2][j] / hx) / (beta_2 + beta_1 / hx)

        u2 = np.zeros((lx, ly))

        for i in range(lx - 1):
            aa, bb, cc, dd = np.zeros(lx), np.zeros(lx), np.zeros(lx), np.zeros(lx)

            bb[0] = hy * beta_3 - alpha_3
            cc[0] = alpha_3
            dd[0] = phi_3() * hy

            aa[-1] = - alpha_4
            bb[-1] = hy * beta_4 + alpha_4
            dd[-1] = phi_4(x[i], t[k]) * hy

            for j in range(1, ly - 1):
                aa[j] = b
                bb[j] = - 2 * (hy ** 2) / tau - 2 * b
                cc[j] = b
                dd[j] = -2 * (hy ** 2) * u1[i][j] / tau \
                        - a * (hy ** 2) * (u1[i+1][j] - 2 * u1[i][j] + u1[i-1][j]) / (hx ** 2) \
                        - (hy ** 2) * f(x[i], y[j], t_half)

            xx = tridiagonal_matrix_algorithm(aa, bb, cc, dd)

            for j in range(ly):
                u2[i][j] = xx[j]
                u2[0][j] = (phi_1() - alpha_1 * u2[1][j] / hx) / (alpha_2 - alpha_1 / hx)
                u2[-1][j] = (phi_2(y[j], t[k]) + beta_1 * u2[-2][j] / hx) / (beta_2 + beta_1 / hx)

        for i in range(lx):
            u2[i][0] = (phi_3() - alpha_3 * u2[i][1] / hy) / (beta_3 - alpha_3 / hy)
            u2[i][-1] = (phi_4(x[i], t[k]) + alpha_4 * u2[i][-2] / hy) / (beta_4 + alpha_4 / hy)

        for i in range(lx):
            for j in range(ly):
                u[i][j][k] = u2[i][j]

    return u


def mds():
    u = np.zeros((lx, ly, lt))

    for i in range(lx):
        for j in range(ly):
            u[i][j][0] = psi()

    for k in range(1, lt):
        u1 = np.zeros((lx, ly))
        t_half = t[k] - tau / 2
        u_k_1 = u[:, :, k - 1]

        for j in range(ly - 1):
            aa, bb, cc, dd = np.zeros(lx), np.zeros(lx), np.zeros(lx), np.zeros(lx)

            bb[0] = hx * alpha_2 - alpha_1
            cc[0] = alpha_1
            dd[0] = phi_1() * hx

            aa[-1] = - beta_1
            bb[-1] = hx * beta_2 + beta_1
            dd[-1] = phi_2(y[j], t_half) * hx

            for i in range(1, lx - 1):
                aa[i] = a
                bb[i] = - (hx ** 2) / tau - 2 * a
                cc[i] = a
                dd[i] = - (hx ** 2) * u_k_1[i][j] / tau - (hx ** 2) * f(x[i], y[j], t_half) / 2

            xx = tridiagonal_matrix_algorithm(aa, bb, cc, dd)

            for i in range(lx):
                u1[i][j] = xx[i]
                u1[i][0] = (phi_3() - alpha_3 * u1[i][1] / hy) / (beta_3 - alpha_3 / hy)
                u1[i][-1] = (phi_4(x[i], t_half) + alpha_4 * u1[i][-2] / hy) / (beta_4 + alpha_4 / hy)

        for j in range(ly):
            u1[0][j] = (phi_1() - alpha_1 * u1[1][j] / hx) / (alpha_2 - alpha_1 / hx)
            u1[-1][j] = (phi_2(y[j], t_half) + beta_1 * u1[-2][j] / hx) / (beta_2 + beta_1 / hx)

        u2 = np.zeros((lx, ly))

        for i in range(lx - 1):
            aa, bb, cc, dd = np.zeros(lx), np.zeros(lx), np.zeros(lx), np.zeros(lx)

            bb[0] = hy * beta_3 - alpha_3
            cc[0] = alpha_3
            dd[0] = phi_3() * hy

            aa[-1] = - alpha_4
            bb[-1] = hy * beta_4 + alpha_4
            dd[-1] = phi_4(x[i], t[k]) * hy

            for j in range(1, ly - 1):
                aa[j] = b
                bb[j] = - (hy ** 2) / tau - 2 * b
                cc[j] = b
                dd[j] = - (hy ** 2) * u1[i][j] / tau - (hy ** 2) * f(x[i], y[j], t[k]) / 2

            xx = tridiagonal_matrix_algorithm(aa, bb, cc, dd)

            for j in range(ly):
                u2[i][j] = xx[j]
                u2[0][j] = (phi_1() - alpha_1 * u2[1][j] / hx) / (alpha_2 - alpha_1 / hx)
                u2[-1][j] = (phi_2(y[j], t[k]) + beta_1 * u2[-2][j] / hx) / (beta_2 + beta_1 / hx)

        for i in range(lx):
            u2[i][0] = (phi_3() - alpha_3 * u2[i][1] / hy) / (beta_3 - alpha_3 / hy)
            u2[i][-1] = (phi_4(x[i], t[k]) + alpha_4 * u2[i][-2] / hy) / (beta_4 + alpha_4 / hy)

        for i in range(lx):
            for j in range(ly):
                u[i][j][k] = u2[i][j]

    return u


def solve_error(u):
    error = []
    u_new = u[:, :, kt]
    for i in range(ly):
        diff = analytical_solution(x, y[i], t[kt]) - u_new[:, i]

        norm = 0
        for j in range(len(diff)):
            norm += diff[j] ** 2

        error = np.append(error, np.sqrt(norm))

    return error / np.sqrt(ly)


def print_parameters():
    print(f"nx: {num_x}, ny: {num_y}")
    print(f"hx: {hx}, hy: {hy}")
    print(f"tau: {tau}")
    print(f"x[kt]: {x[kt]}, y[kt]: {y[kt]}, t[kt]: {t[kt]}")
