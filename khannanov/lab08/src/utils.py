import numpy as np
import math
import matplotlib.pyplot as plt
import random
from schema import *

class Data:
	def __init__(self, params):
		self.a = params['a']
		self.b = params['b']
		self.c = params['c']
		self.d = params['d']
		self.lx = params['lx']
		self.ly = params['ly']
		self.f = params['f']
		self.alpha1 = params['alpha1']
		self.alpha2 = params['alpha2']
		self.beta1 = params['beta1']
		self.beta2 = params['beta2']
		self.gamma1 = params['gamma1']
		self.gamma2 = params['gamma2']
		self.delta1 = params['delta1']
		self.delta2 = params['delta2']
		self.phi11 = params['phi11']
		self.phi21 = params['phi21']
		self.phi12 = params['phi12']
		self.phi22 = params['phi22']
		self.psi = params['psi']
		self.solution = params['solution']
            

def psi_0(x, t):
    return 0.0

def psi_1(x, t):
    return x * math.cos(t)

def phi_0(y, t):
    return 0.0

def phi_1(y, t):
    return y * math.cos(t)

def u0(x, y):
    return x*y

def u(x, y, t):
    return x*y * math.cos(t)
		

def solve_tridiagonal_system(A, B, C, D):
    size = len(A)
    P, Q = [], []
    P.append(-C[0] / B[0])
    Q.append(D[0] / B[0])

    for i in range(1, size):
        P_tmp = -C[i] / (B[i] + A[i] * P[i - 1])
        Q_tmp = (D[i] - A[i] * Q[i - 1]) / (B[i] + A[i] * P[i - 1])
        P.append(P_tmp)
        Q.append(Q_tmp)

    X = [0 for _ in range(size)]
    X[size - 1] = Q[size - 1]

    for i in range(size - 2, -1, -1):
        X[i] = P[i] * X[i + 1] + Q[i]

    return X

def max_difference(L, u, nx, ny):
    mx_difference = 0
    for i in range(nx):
        for j in range(ny):
            mx_difference = max(mx_difference, abs(u[i][j] - L[i][j]))
    return mx_difference


def compareError(a, b):
	err = 0
	lst = [abs(i - j) for i, j in zip(a, b)]
	for each in lst:
		err = max(err, each)
	return err


def approximation(x_list, y, t, dict_):
	res = []
	for i in range(len(x_list)):
		res.append(dict_['numerical'][i][y][t])
	return res


def solving(x_list, y, t):
	res = []
	for xi in x_list:
		res.append(xi * y * np.cos(t))
	return res

def RealZByTime(lx0, lx1, ly0, ly1, t, f):
    x = np.arange(lx0, lx1 + 0.002, 0.002)
    y = np.arange(ly0, ly1 + 0.002, 0.002)
    X = np.ones((y.shape[0], x.shape[0]))
    Y = np.ones((x.shape[0], y.shape[0]))
    Z = np.ones((y.shape[0], x.shape[0]))
    for i in range(Y.shape[0]):
        Y[i] = y
    Y = Y.T
    for i in range(X.shape[0]):
        X[i] = x
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            Z[i, j] = f(X[i, j], Y[i, j], t)
    return X, Y, Z

def Error(X, Y, t, z, ut = u):
    ans = 0.0
    for i in range(len(z)):
        for j in  range(len(z[i])):
            ans = max(abs(ut(X[i][j], Y[i][j], t) - z[i][j]), ans)
    return (ans / len(z) / len(z[0]))


def StepSlice(lst, step):
    return lst[step]


def PlotByTime(X, Y, T, Z, j, extrems, plot_true=True):
	t = T[j]
	z = Z[j]
	fig = plt.figure(num=1, figsize=(20, 13), clear=True)
	ax = fig.add_subplot(1, 1, 1, projection='3d')
	ax.plot_surface(np.array(X), np.array(Y), np.array(z))
	if plot_true:
		ax.plot_wireframe(*RealZByTime(0, 1, 0, 1, t, u), color="green")
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	ax.set_title(
		't = ' + str(round(t, 8)) + " error = " + str(round(Error(X, Y, t, z), 11)),
		loc="right", fontsize=25
	)
	ax.set_zlim(extrems[0], extrems[1])
	fig.tight_layout()
	plt.close(fig)
	return fig


def SquareMinMax(z):
	minimum, maximum = z[0][0], z[0][0]
	for i in range(len(z)):
		for j in range(len(z[i])):
			minimum = z[i][j] if z[i][j] < minimum else minimum
			maximum = z[i][j] if z[i][j] > maximum else maximum
	return minimum, maximum


def SearchMinMax(zz):
	minimum, maximum = 0.0, 0.0
	for z in zz:
		minmax = SquareMinMax(z)
		minimum = minmax[0] if minmax[0] < minimum else minimum
		maximum = minmax[1] if minmax[1] > maximum else maximum
	return minimum, maximum


def GetGraphicH(solver, time = 0, tsteps = 40):
    h, e = [], []
    for N in range(4, 20, 1):
        x, y, t, z = solver(Nx = N, Ny = N, K = tsteps)
        h.append(solver.hx)
        e.append(Error(x, y, t[time], z[time]))
    return h, e

first = Schema(T = 2*math.pi, order2nd = False)
second = Schema(T = 2*math.pi, order2nd = True)

TSTEPS = 100
time = random.randint(0, TSTEPS - 1)
h1, e1 = GetGraphicH(first, time, TSTEPS)
h2, e2 = GetGraphicH(second, time, TSTEPS)

def GetGraphicTau(solver):
    tau = []
    e = []
    for K in range(15, 100, 2):
        x, y, t, z = solver(Nx = 10, Ny = 10, K = K)
        tau.append(solver.tau)
        time = K // 2
        e.append(Error(x, y, t[time], z[time]))
    return tau, e

tau1, e1 = GetGraphicTau(first)
tau2, e2 = GetGraphicTau(second)

def FullError(X, Y, T, Z):
    ans = 0.0
    for k in range(len(T)):
        for i in range(len(X)):
            for j in range(len(X[i])):
                ans = max(abs(u(X[i][j], Y[i][j], T[k]) - Z[k][i][j]), ans)
    return (ans / len(T) / len(X) / len(X[0]))


TimeList = [random.randint(0, 40 - 1) for i in range(4)]