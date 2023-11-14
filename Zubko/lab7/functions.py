import numpy as np


def analytic_solution(x, y):
    return np.exp(-x) * np.cos(x) * np.cos(y)


def phi1(y):
    return np.cos(y)


def phi2(y):
    return 0


def phi3(x):
    return np.exp(-x) * np.cos(x)


def phi4(x):
    return 0
