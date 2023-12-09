import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def tma(a, b, c, d):
    size = len(a)
    p, q = [], []
    p.append(-c[0] / b[0])
    q.append(d[0] / b[0])

    for i in range(1, size):
        p_tmp = -c[i] / (b[i] + a[i] * p[i - 1])
        q_tmp = (d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1])
        p.append(p_tmp)
        q.append(q_tmp)

    x = [0 for _ in range(size)]
    x[size - 1] = q[size - 1]

    for i in range(size - 2, -1, -1):
        x[i] = p[i] * x[i + 1] + q[i]

    return x


def norm_inf(A):
    n = len(A)
    norm = 0
    for i in range(n):
        sum_ = 0
        for j in range(n):
            sum_ += abs(A[i][j])
        norm = sum_ if norm < sum_ else norm
    return norm


def norm_inf_vec(A):
    n = len(A)
    norm = 0
    for i in range(n):
        if abs(A[i]) > norm:
            norm = abs(A[i])
    return norm

def diff(L, u, nx, ny):
	mx = 0
	for i in range(nx):
		for j in range(ny):
			mx = max(mx, abs(u[i][j] - L[i][j]))
	return mx


class EquationParameters:
    def __init__(self, parameters):
        for key, value in parameters.items():
            setattr(self, key, value)


class ParabolicSolver:
	def __init__(self, params, equation_type):
		self.data = EquationParameters(params)
		self.hx = 0
		self.hy = 0
		self.tau = 0
		try:
			self.solve_func = getattr(self, f'{equation_type}_solver')
		except:
			raise Exception("This type does not exist")

	def getCoeffs(self, n):
		aa = np.zeros(len(n))
		bb = np.zeros(len(n))
		cc = np.zeros(len(n))
		dd = np.zeros(len(n))

		return aa, bb, cc, dd

	def computeCoeffs(self, x, y, t2, j):
		aa, bb, cc, dd = self.getCoeffs(x)
		bb[0] = self.hx * self.data.alpha2 - self.data.alpha1
		bb[-1] = self.hx * self.data.beta2 + self.data.beta1
		cc[0] = self.data.alpha1
		aa[-1] = -self.data.beta1
		dd[0] = self.data.phi11(y[j], t2) * self.hx
		dd[-1] = self.data.phi12(y[j], t2) * self.hx

		return aa, bb, cc, dd

	def prepare(self, nx, ny, T, K):
		self.hx = self.data.lx / nx
		self.hy = self.data.ly / ny
		self.tau = T / K
		x = np.arange(0, self.data.lx + self.hx, self.hx)
		y = np.arange(0, self.data.ly + self.hy, self.hy)
		t = np.arange(0, T + self.tau, self.tau)

		return x, y, t

	def initalizeU(self, x, y, t):
		u = np.zeros((len(x), len(y), len(t)))
		for i in range(len(x)):
			for j in range(len(y)):
				u[i][j][0] = self.data.psi(x[i], y[j])

		return u

	def solve(self, nx, ny, T, K):
		self.hx = self.data.lx / nx
		self.hy = self.data.ly / ny
		self.tau = T / K

		x, y, t = self.prepare(nx, ny, T, K)

		uu = self.initalizeU(x, y, t)

		return self.solve_func(x, y, t, uu)

	def analyticSolve(self, nx, ny, T, K):
		x, y, t = self.prepare(nx, ny, T, K)

		uu = np.zeros((len(x), len(y), len(t)))

		for i in range(len(x)):
			for j in range(len(y)):
				for k in range(len(t)):
					uu[i][j][k] = self.data.solution(x[i], y[j], t[k])

		return uu

	def parallelDirections_solver(self, x, y, t, uu):
		for k in range(1, len(t)):
			u1 = np.zeros((len(x), len(y)))
			t2 = t[k] - self.tau / 2
			for j in range(len(y) - 1):
				aa, bb, cc, dd = self.computeCoeffs(x, y, t2, j)
				for i in range(len(x) - 1):
					aa[i] = self.data.a - self.hx * self.data.c / 2
					bb[i] = self.hx ** 2 - 2 * (self.hx ** 2) / self.tau - 2 * self.data.a
					cc[i] = self.data.a + self.hx * self.data.c / 2
					dd[i] = -2 * (self.hx ** 2) * uu[i][j][k - 1] / self.tau
					- self.data.b * (self.hx ** 2) * (uu[i][j + 1][k - 1]
													  - 2 * uu[i][j][k - 1] + uu[i][j - 1][k - 1]) / (self.hy ** 2)
					- self.data.d * (self.hx ** 2) * (uu[i][j + 1][k - 1] - uu[i][j - 1][k - 1]) / (2 * self.hy ** 2)
					- (self.hx ** 2) * self.data.f(x[i], y[j], t[k])

				xx = tma(aa, bb, cc, dd)
				for i in range(len(x)):
					u1[i][j] = xx[i]
					u1[i][0] = (self.data.phi21(x[i], t2) - self.data.gamma1 * u1[i][1] / self.hy) / (
							self.data.gamma2 - self.data.gamma1 / self.hy)
					u1[i][-1] = (self.data.phi22(x[i], t2) + self.data.delta1 * u1[i][-2] / self.hy) / (
							self.data.delta2 + self.data.delta1 / self.hy)
			for j in range(len(y)):
				u1[0][j] = (self.data.phi11(y[j], t2) - self.data.alpha1 * u1[1][j] / self.hx) / (
							self.data.alpha2 - self.data.alpha1 / self.hx)
				u1[-1][j] = (self.data.phi12(y[j], t2) + self.data.beta1 * u1[-2][j] / self.hx) / (
							self.data.beta2 + self.data.beta1 / self.hx)
			####
			u2 = np.zeros((len(x), len(y)))
			for i in range(len(x) - 1):
				aa, bb, cc, dd = self.getCoeffs(y)
				bb[0] = self.hy * self.data.gamma2 - self.data.gamma1
				bb[-1] = self.hy * self.data.delta2 + self.data.delta1
				cc[0] = self.data.gamma1
				aa[-1] = -self.data.delta1
				dd[0] = self.data.phi21(x[i], t[k]) * self.hy
				dd[-1] = self.data.phi22(x[i], t[k]) * self.hy

				for j in range(len(y) - 1):
					aa[j] = self.data.b - self.hy * self.data.d / 2
					bb[j] = self.hy ** 2 - 2 * (self.hy ** 2) / self.tau - 2 * self.data.b
					cc[j] = self.data.b + self.hy * self.data.d / 2
					dd[j] = -2 * (self.hy ** 2) * u1[i][j] / self.tau
					- self.data.a * (self.hy ** 2) * (u1[i + 1][j]
													  - 2 * u1[i][j] + u1[i - 1][j]) / (self.hx ** 2)
					- self.data.c * (self.hy ** 2) * (u1[i + 1][j] - u1[i - 1][j]) / (2 * self.hx ** 2)
					- (self.hy ** 2) * self.data.f(x[i], y[j], t[k])
				xx = tma(aa, bb, cc, dd)
				for j in range(len(y)):
					u2[i][j] = xx[j]
					u2[0][j] = (self.data.phi11(y[j], t[k]) - self.data.alpha1 * u2[1][j] / self.hx) / (
								self.data.alpha2 - self.data.alpha1 / self.hx)
					u2[-1][j] = (self.data.phi12(y[j], t[k]) + self.data.beta1 * u2[-2][j] / self.hx) / (
								self.data.beta2 + self.data.beta1 / self.hx)
			for i in range(len(x)):
				u2[i][0] = (self.data.phi21(x[i], t[k]) - self.data.gamma1 * u2[i][1] / self.hy) / (
							self.data.gamma2 - self.data.gamma1 / self.hy)
				u2[i][-1] = (self.data.phi22(x[i], t[k]) + self.data.delta1 * u2[i][-2] / self.hy) / (
							self.data.delta2 + self.data.delta1 / self.hy)
			for i in range(len(x)):
				for j in range(len(y)):
					uu[i][j][k] = u2[i][j]
		return uu

	def fractionalSteps_solver(self, x, y, t, uu):
		for k in range(len(t)):
			u1 = np.zeros((len(x), len(y)))
			t2 = t[k] - self.tau / 2
			for j in range(len(y) - 1):
				aa, bb, cc, dd = self.computeCoeffs(x, y, t2, j)
				for i in range(len(x) - 1):
					aa[i] = self.data.a
					bb[i] = -(self.hx ** 2) / self.tau - 2 * self.data.a
					cc[i] = self.data.a
					dd[i] = -(self.hx ** 2) * uu[i][j][k - 1] / self.tau - (self.hx ** 2) * self.data.f(x[i], y[j],
																										t2) / 2
				xx = tma(aa, bb, cc, dd)
				for i in range(len(x)):
					u1[i][j] = xx[i]
					u1[i][0] = (self.data.phi21(x[i], t2) - self.data.gamma1 * u1[i][1] / self.hy) / (
							self.data.gamma2 - self.data.gamma1 / self.hy)
					u1[i][-1] = (self.data.phi22(x[i], t2) + self.data.delta1 * u1[i][-2] / self.hy) / (
							self.data.delta2 + self.data.delta1 / self.hy)
			for j in range(len(y)):
				u1[0][j] = (self.data.phi11(y[j], t2) - self.data.alpha1 * u1[1][j] / self.hx) / (
						self.data.alpha2 - self.data.alpha1 / self.hx)
				u1[-1][j] = (self.data.phi12(y[j], t2) + self.data.beta1 * u1[-2][j] / self.hx) / (
						self.data.beta2 + self.data.beta1 / self.hx)
			#####
			u2 = np.zeros((len(x), len(y)))
			for i in range(len(x) - 1):
				aa, bb, cc, dd = self.getCoeffs(y)
				bb[0] = self.hy * self.data.gamma2 - self.data.gamma1
				bb[-1] = self.hy * self.data.delta2 + self.data.delta1
				cc[0] = self.data.gamma1
				aa[-1] = -self.data.delta1
				dd[0] = self.data.phi21(x[i], t[k]) * self.hy
				dd[-1] = self.data.phi22(x[i], t[k]) * self.hy

				for j in range(len(y) - 1):
					aa[j] = self.data.b
					bb[j] = -(self.hy ** 2) / self.tau - 2 * self.data.b
					cc[j] = self.data.b
					dd[j] = -(self.hy ** 2) * u1[i][j] / self.tau - (self.hy ** 2) * self.data.f(x[i], y[j], t[k]) / 2
				xx = tma(aa, bb, cc, dd)
				for j in range(len(y)):
					u2[i][j] = xx[j]
					u2[0][j] = (self.data.phi11(y[j], t[k]) - self.data.alpha1 * u2[1][j] / self.hx) / (
							self.data.alpha2 - self.data.alpha1 / self.hx)
					u2[-1][j] = (self.data.phi12(y[j], t[k]) + self.data.beta1 * u2[-2][j] / self.hx) / (
							self.data.beta2 + self.data.beta1 / self.hx)
			for i in range(len(x)):
				u2[i][0] = (self.data.phi21(x[i], t[k]) - self.data.gamma1 * u2[i][1] / self.hy) / (
						self.data.gamma2 - self.data.gamma1 / self.hy)
				u2[i][-1] = (self.data.phi22(x[i], t[k]) + self.data.delta1 * u2[i][-2] / self.hy) / (
						self.data.delta2 + self.data.delta1 / self.hy)
			for i in range(len(x)):
				for j in range(len(y)):
					uu[i][j][k] = u2[i][j]
		return uu


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


def draw_u_x(dict_, nx, ny, T, K, time, save_file="plot_u_x.png"):
	fig = plt.figure()
	hx = 1 / nx
	hy = 1 / ny
	tau = T / K

	x = np.arange(0, 1 + hx, hx)
	y = np.arange(0, 1 + hy, hy)
	t = np.arange(0, T + tau, tau)

	#     z1 = np.array(dict_['numerical'])
	#     z2 = np.array(dict_['analytic'])

	z1 = []
	for i in range(len(x)):
		z1.append(approximation(x, 10, time, dict_))

	z2 = solving(x, y[10], t[time])

	plt.title('U from x')
	plt.plot(x, z1[time], color='r', label='numerical')
	#     print(x)
	#     print(z1)
	plt.plot(x, z2, color='b', label='analytic')
	plt.legend(loc='best')
	plt.ylabel('U')
	plt.xlabel('x')
	plt.savefig(save_file)
	plt.show()

	err = []
	print(f"len={len(z1[time])}\n{z1[time]}")
	print()
	print(f"len={len(z2)}\n{z2}")
	for i in range(len(z1[time])):
		err.append(compareError(z1[time], z2))
	plt.title('Error from y')
	plt.plot(y, err, color='b', label='err')
	plt.legend(loc='best')
	plt.ylabel('Err')
	plt.xlabel('t')
	plt.savefig('err.png')
	plt.show()


data = {'equation_type': 'parallelDirections', 'nx': 40, 'ny': 40, 'T': 5, 'K': 200}

if __name__ == '__main__':
	equation_type = data['equation_type']
	nx, ny, T, K = int(data['nx']), int(data['ny']), int(data['T']), int(data['K'])
	params = {
		'a': 1,
		'b': 1,
		'c': 0,
		'd': 0,
		'lx': 1,
		'ly': 1,
		'f': lambda x, y, t: - x * y * np.sin(t),
		'alpha1': 0,
		'alpha2': 1,
		'beta1': 0,
		'beta2': 1,
		'gamma1': 0,
		'gamma2': 1,
		'delta1': 0,
		'delta2': 1,
		'phi11': lambda y, t: 0,
		'phi12': lambda y, t: y * np.cos(t),
		'phi21': lambda x, t: 0,
		'phi22': lambda x, t: x * np.cos(t),
		'psi': lambda x, y: x * y,
		'solution': lambda x, y, t: x * y * np.cos(t),
	}

	solver = ParabolicSolver(params, equation_type)
	solved = solver.solve(nx, ny, T, K)
	ans = {
		'numerical': solved,
		'analytic': solver.analyticSolve(nx, ny, T, K)
	}

	draw_u_x(ans, nx, ny, T, K, 30)