import numpy as np

from task import a, c


def analytic_solution(x, t):
    return np.exp((c - a) * t) * np.sin(x)


def psi(x):
    return np.sin(x)


def phi0(t):
    return np.exp((c - a) * t)


def phi1(t):
    return np.exp((c - a) * t)
