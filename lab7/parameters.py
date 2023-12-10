import numpy as np

bx = 2
by = 2
c = 4

alpha_1 = 0
beta_1 = 1

alpha_2 = 0
beta_2 = 1

alpha_3 = 0
beta_3 = 1

alpha_4 = 0
beta_4 = 1

lx = np.pi / 2
ly = np.pi / 2

K = 50
N = 50
X = np.linspace(0, lx, K)

hx = lx / K
hy = ly / N
U = []

tt = 5

omega = 1.9
eps = 0.00001


def phi_1(y):
    return np.exp(-y) * np.cos(y)


def phi_2(y):
    return 0


def phi_3(x):
    return np.exp(-x) * np.cos(x)


def phi_4(x):
    return 0


def analytical_solution(x, y):
    return np.exp(-x - y) * np.cos(x) * np.cos(y)