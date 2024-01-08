import math


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

# analytic solve
def u(x, y, t):
    return x*y * math.cos(t)


class Schema:
	def __init__(self, rho=u0, psi0=psi_0, psi1=psi_1, phi0=phi_0, phi1=phi_1,
				 lx0=0, lx1=1.0, ly0=0, ly1=1.0, T=3, order2nd=True):
		self.psi0 = psi0
		self.psi1 = psi1
		self.phi0 = phi0
		self.phi1 = phi1
		self.rho0 = rho
		self.T = T
		self.lx0 = lx0
		self.lx1 = lx1
		self.ly0 = ly0
		self.ly1 = ly1
		self.tau = None
		self.hx = None
		self.hy = None
		self.order = order2nd
		self.Nx = None
		self.Ny = None
		self.K = None
		self.cx = None
		self.bx = None
		self.cy = None
		self.by = None
		self.hx2 = None
		self.hy2 = None

	def set_l0_l1(self, lx0, lx1, ly0, ly1):
		self.lx0 = lx0
		self.lx1 = lx1
		self.ly0 = ly0
		self.ly1 = ly1

	def set_T(self, T):
		self.T = T

	def CalculateH(self):
		self.hx = (self.lx1 - self.lx0) / (self.Nx - 1)
		self.hy = (self.ly1 - self.ly0) / (self.Ny - 1)
		self.hx2 = self.hx * self.hx
		self.hy2 = self.hy * self.hy

	def CalculateTau(self):
		self.tau = self.T / (self.K - 1)

	@staticmethod
	def race_method(A, b):
		P = [-item[2] for item in A]
		Q = [item for item in b]
		P[0] /= A[0][1]
		Q[0] /= A[0][1]
		for i in range(1, len(b)):
			z = (A[i][1] + A[i][0] * P[i - 1])
			P[i] /= z
			Q[i] -= A[i][0] * Q[i - 1]
			Q[i] /= z
		for i in range(len(Q) - 2, -1, -1):
			Q[i] += P[i] * Q[i + 1]
		return Q

	@staticmethod
	def nparange(start, end, step=1):
		now = start
		e = 0.00000000001
		while now - e <= end:
			yield now
			now += step

	def CalculateLeftEdge(self, X, Y, t, square):
		for i in range(self.Ny):
			square[i][0] = self.phi0(Y[i][0], t)

	def CalculateRightEdge(self, X, Y, t, square):
		for i in range(self.Ny):
			square[i][-1] = self.phi1(Y[i][-1], t)

	def CalculateBottomEdge(self, X, Y, t, square):
		for j in range(1, self.Nx - 1):
			square[0][j] = self.psi0(X[0][j], t)

	def CalculateTopEdge(self, X, Y, t, square):
		for j in range(1, self.Nx - 1):
			square[-1][j] = self.psi1(X[-1][j], t)

	def CalculateLineFirstStep(self, i, X, Y, t, last_square, now_square):
		hy2 = self.hy2
		hx2 = self.hx2
		b = self.bx
		c = self.cx
		A = [(0, b, c)]
		w = [
			-self.cy * self.order * last_square[i - 1][1] -
			((self.order + 1) * hx2 * hy2 - 2 * self.cy * self.order) * last_square[i][1] -
			self.cy * self.order * last_square[i + 1][1] +
			self.tau * hy2 * hx2 * X[i][1] * Y[i][1] * math.sin(t) -
			c * now_square[i][0]
		]
		A.extend([(c, b, c) for _ in range(2, self.Nx - 2)])
		w.extend([
			-self.cy * self.order * last_square[i - 1][j] -
			((self.order + 1) * hx2 * hy2 - 2 * self.cy * self.order) * last_square[i][j] -
			self.cy * self.order * last_square[i + 1][j] +
			self.tau * hy2 * hx2 * X[i][j] * Y[i][j] * math.sin(t)
			for j in range(2, self.Nx - 2)
		])
		A.append((c, b, 0))
		w.append(
			-self.cy * self.order * last_square[i - 1][-2] -
			((self.order + 1) * hx2 * hy2 - 2 * self.cy * self.order) * last_square[i][-2] -
			self.cy * self.order * last_square[i + 1][-2] +
			self.tau * hy2 * hx2 * X[i][-2] * Y[i][-2] * math.sin(t) -
			c * now_square[i][-1]
		)
		line = self.race_method(A, w)
		for j in range(1, self.Nx - 1):
			now_square[i][j] = line[j - 1]

	def CalculateLineSecondStep(self, j, X, Y, t, last_square, now_square):
		hx2 = self.hx2
		hy2 = self.hy2
		c = self.cy
		b = self.by
		A = [(0, b, c)]
		w = [
			-self.cx * self.order * last_square[1][j - 1] -
			((self.order + 1) * hx2 * hy2 - 2 * self.cx * self.order) * last_square[1][j] -
			self.cx * self.order * last_square[1][j + 1] +
			self.tau * hy2 * hx2 * X[1][j] * Y[1][j] * math.sin(t) -
			c * now_square[0][j]
		]
		A.extend([(c, b, c) for _ in range(2, self.Ny - 2)])
		w.extend([
			-self.cx * self.order * last_square[i][j - 1] -
			((self.order + 1) * hx2 * hy2 - 2 * self.cx * self.order) * last_square[i][j] -
			self.cx * self.order * last_square[i][j + 1] +
			self.tau * hy2 * hx2 * X[i][j] * Y[i][j] * math.sin(t)
			for i in range(2, self.Ny - 2)
		])
		A.append((c, b, 0))
		w.append(
			-self.cx * self.order * last_square[-2][j - 1] -
			((self.order + 1) * hx2 * hy2 - 2 * self.cx * self.order) * last_square[-2][j] -
			self.cx * self.order * last_square[-2][j + 1] +
			self.tau * hy2 * hx2 * X[-2][j] * Y[-2][j] * math.sin(t) -
			c * now_square[-1][j]
		)
		line = self.race_method(A, w)
		for i in range(1, self.Ny - 1):
			now_square[i][j] = line[i - 1]

	def CalculateSquare(self, X, Y, t, last_square):
		square = [[0.0 for _ in range(self.Nx)] for _ in range(self.Ny)]
		self.CalculateLeftEdge(X, Y, t - 0.5 * self.tau, square)
		self.CalculateRightEdge(X, Y, t - 0.5 * self.tau, square)
		self.CalculateBottomEdge(X, Y, t - 0.5 * self.tau, square)
		self.CalculateTopEdge(X, Y, t - 0.5 * self.tau, square)
		for i in range(1, self.Ny - 1):
			self.CalculateLineFirstStep(i, X, Y, t - 0.5 * self.tau, last_square, square)
		last_square = square
		square = [[0.0 for _ in range(self.Nx)] for _ in range(self.Ny)]
		self.CalculateLeftEdge(X, Y, t, square)
		self.CalculateRightEdge(X, Y, t, square)
		self.CalculateBottomEdge(X, Y, t, square)
		self.CalculateTopEdge(X, Y, t, square)
		for j in range(1, self.Nx - 1):
			self.CalculateLineSecondStep(j, X, Y, t, last_square, square)
		return square

	def init_t0(self, X, Y):
		first = [[0.0 for _ in range(self.Nx)] for _ in range(self.Ny)]
		for i in range(self.Ny):
			for j in range(self.Nx):
				first[i][j] = self.rho0(X[i][j], Y[i][j])
		return first

	def __call__(self, Nx=20, Ny=20, K=20):
		self.Nx, self.Ny, self.K = Nx, Ny, K
		self.CalculateTau()
		self.CalculateH()

		self.bx = -2 * self.tau * self.hy2
		self.bx -= (1 + self.order) * self.hx2 * self.hy2
		self.cx = self.tau * self.hy2

		self.cy = self.tau * self.hx2
		self.by = -2 * self.tau * self.hx2
		self.by -= (1 + self.order) * self.hx2 * self.hy2
		x = list(self.nparange(self.lx0, self.lx1, self.hx))
		y = list(self.nparange(self.ly0, self.ly1, self.hy))
		X = [x for _ in range(self.Ny)]
		Y = [[y[i] for _ in x] for i in range(self.Ny)]

		taus = [0.0]
		ans = [self.init_t0(X, Y)]
		for t in self.nparange(self.tau, self.T, self.tau):
			ans.append(self.CalculateSquare(X, Y, t, ans[-1]))
			taus.append(t)
		return X, Y, taus, ans