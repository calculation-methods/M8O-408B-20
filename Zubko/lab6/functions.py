import numpy as np


def analytic_solution(x, t):
    return np.exp(-t - x) * np.sin(x) * np.sin(2 * t)


def phi0(t):
    return 0


def phi1(t):
    return 0


def psi1(x):
    return 0


def psi2(x):
    return 2 * np.exp(-x) * np.sin(x)
