import numpy as np

kx = 5
kt = 5

# 1 ---
a = 1
b = 1
mu = 1
# -----

# # 2 ---
# a = 2
# b = 1
# mu = 1
# # -----

# # 3 ---
# a = 1
# b = 2
# mu = 1
# # -----

# # 4 ---
# a = 1
# b = 1
# mu = 2
# # -----

alpha_1 = 0
beta_1 = 1

alpha_2 = 1
beta_2 = 0

alpha_3 = 0
beta_3 = 1

alpha_4 = 1
beta_4 = 0

left_border_x = 0
# right_border_x = np.pi
right_border_x = 3.14
num_x = 20
hx = (right_border_x - left_border_x) / num_x

left_border_y = 0
# right_border_y = np.pi
right_border_y = 3.14
num_y = 20
hy = (right_border_y - left_border_y) / num_y

T = 1
K = 50
tau = T / K

x = np.arange(left_border_x, right_border_x + hx, hx)
y = np.arange(left_border_y, right_border_y + hy, hy)
t = np.arange(0, T + tau, tau)

lx, ly, lt = len(x), len(y), len(t)


def analytical_solution(_x, _y, _t):
    return np.sin(_x) * np.sin(_y) * np.sin(mu * _t)


def phi_1():
    return 0


def phi_2(_y, _t):
    return - np.sin(_y) * np.sin(mu * _t)


def phi_3():
    return 0


def phi_4(_x, _t):
    return - np.sin(_x) * np.sin(mu * _t)


def psi():
    return 0


def f(_x, _y, _t):
    return np.sin(_x) * np.sin(_y) * (mu * np.cos(mu * _t) + (a + b) * np.sin(mu * _t))
