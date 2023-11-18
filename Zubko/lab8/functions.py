import numpy as np


def f(x, y, t):
    return -x * y * np.sin(t)


def analytic_solution(x, y, t):
    return x * y * np.cos(t)


def phi1(y, t):
    return 0


def phi2(y, t):
    return 0


def phi3(x, t):
    return 0


def phi4(x, t):
    return 0


def psi(x, y):
    return x * y
