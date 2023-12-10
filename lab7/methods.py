from parameters import *
from copy import deepcopy as dp


def interpolation(u):
    for i in range(1, K - 1):
        for j in range(1, N - 1):
            u[i][j] = u[0][j] + (X[i] - X[0]) * (u[-1][j] - u[0][j]) / (X[-1] - X[0])
    return u


def liebmann_method():
    u = np.zeros((K, N))

    for j in range(0, N):
        u[0][j] = phi_1(hy * j) / beta_1
        u[K-1][j] = phi_2(hy * j) / beta_2

    for i in range(0, K):
        u[i][0] = phi_3(hx * i) / beta_3
        u[i][N-1] = phi_4(hx * i) / beta_4

    u = interpolation(u.copy())

    cnt = 0
    while True:
        cnt += 1
        u_old = dp(u)
        for i in range(1, K - 1):
            for j in range(1, N - 1):

                u[i][j] = (u_old[i+1][j] * (1/hx**2 + bx/(2*hx)) +
                           u_old[i-1][j] * (1/hx**2 - bx/(2*hx)) +
                           u_old[i][j+1] * (1/hy**2 + by/(2*hy)) +
                           u_old[i][j-1] * (1/hy**2 - by/(2*hy))) / (2/hx**2 + 2/hy**2 - c)

        if stop_condition(u_old, u) <= eps:
            break

    return u, cnt


def seidel_method():
    u = np.zeros((K, N))

    for j in range(0, N):
        u[0][j] = phi_1(hy * j) / beta_1
        u[K-1][j] = phi_2(hy * j) / beta_2

    for i in range(0, K):
        u[i][0] = phi_3(hx * i) / beta_3
        u[i][N-1] = phi_4(hx * i) / beta_4

    u = interpolation(u.copy())

    cnt = 0
    while True:
        cnt += 1

        u_old = dp(u)
        for i in range(1, K - 1):
            for j in range(1, N - 1):
                u[i][j] = (u_old[i+1][j] * (1 / hx ** 2 + bx / (2 * hx)) +
                           u[i-1][j] * (1 / hx ** 2 - bx / (2 * hx)) +
                           u_old[i][j+1] * (1 / hy ** 2 + by / (2 * hy)) +
                           u[i][j-1] * (1 / hy ** 2 - by / (2 * hy))) / (2 / hx ** 2 + 2 / hy ** 2 - c)

        if stop_condition(u_old, u) <= eps:
            break

    return u, cnt


def upper_relaxation():
    u = np.zeros((K, N))

    for j in range(0, N):
        u[0][j] = phi_1(hy * j) / beta_1
        u[K-1][j] = phi_2(hy * j) / beta_2

    for i in range(0, K):
        u[i][0] = phi_3(hx * i) / beta_3
        u[i][N-1] = phi_4(hx * i) / beta_4

    u = interpolation(u.copy())

    cnt = 0
    while True:
        cnt += 1
        u_old = dp(u)
        for i in range(1, K - 1):
            for j in range(1, N - 1):
                u[i][j] = u[i][j] + omega * ((u_old[i + 1][j] * (1 / hx ** 2 + bx / (2 * hx)) +
                          u[i - 1][j] * (1 / hx ** 2 - bx / (2 * hx)) +
                          u_old[i][j + 1] * (1 / hy ** 2 + by / (2 * hy)) +
                          u[i][j - 1] * (1 / hy ** 2 - by / (2 * hy))) / (2 / hx ** 2 + 2 / hy ** 2 - c) - u[i][j])

        if stop_condition(u_old, u) <= eps:
            break

    return u, cnt


def stop_condition(u_new, u_old):
    abs_diff = 0
    for i in range(N):
        for j in range(K):
            if abs(u_old[i][j] - u_new[i][j]) > abs_diff:
                abs_diff = abs(u_old[i][j] - u_new[i][j])
    return abs_diff


def solve_error(u):
    error = []

    diff = abs(analytical_solution(X, hy * tt) - u[:, tt])
    error = np.append(error, np.sqrt(diff)) / np.sqrt(len(u))

    return error


def print_parameters():
    print(f"K: {K}, N: {N}")
    print(f"hx: {hx}, hy: {hy}")
    print(f"x[tt]: {X[tt]}")

