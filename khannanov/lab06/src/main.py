import numpy as np
import matplotlib.pyplot as plt


class HyperbolicSolver:
    def __init__(self, parameters):
        for name, value in parameters.items():
            setattr(self, name, value)
        first = parameters['algorithm'][0].upper()
        word = parameters['algorithm'][1:].lower()
        function_name = first + word

    def solve(self, N, K, T):
        self.h = self.l / N
        self.tau = T / K
        self.sigma = (self.tau ** 2) / (self.h ** 2)
        return self.function_name(N, K, T)

    def analytical_solution_matrix(self, N, K, T):
        self.h = self.l / N
        self.tau = T / K
        self.sigma = (self.tau ** 2) / (self.h ** 2)
        self.u = np.zeros((K, N))
        for k in range(K):
            for j in range(N):
                self.u[k][j] = self.analytical_solution(j * self.h, k * self.tau)
        return self.u

    def run_through_method(self):
        size = len(self.a)
        p = np.zeros(size)
        q = np.zeros(size)
        p[0] = (-self.c[0] / self.b[0])
        q[0] = (self.d[0] / self.b[0])

        for i in range(1, size):
            p[i] = -self.c[i] / (self.b[i] + self.a[i] * p[i - 1])
            q[i] = (self.d[i] - self.a[i] * q[i - 1]) / (self.b[i] + self.a[i] * p[i - 1])

        x = np.zeros(size)
        x[-1] = q[-1]

        for i in range(size - 2, -1, -1):
            x[i] = p[i] * x[i + 1] + q[i]

        return x

    def Implicit(self, N, K, T):
        self.u = np.zeros((K, N))
        for j in range(N):
            x_current = j * self.h
            self.u[0][j] = self.psi1(x_current)
            if self.accuracy_level == 1:
                self.u[1][j] = self.psi1(x_current) + self.psi2(x_current) * self.tau + self.psi12(x_current) * self.tau ** 2 / 2
            elif self.accuracy_level == 2:
                k = self.tau ** 2 / 2
                self.u[1][j] = (1 + self.c_const * k) * self.psi2(x_current) + self.a_const * k * self.psi12(x_current) + self.b_const * k * self.psi11(x_current) + (self.tau - self.d_const * k) * self.psi1(x_current) + k * self.function()

        self.a = np.zeros(N)
        self.b = np.zeros(N)
        self.c = np.zeros(N)
        self.d = np.zeros(N)
        for k in range(2, K):
            for j in range(1, N - 1):
                self.a[j] = self.sigma
                self.b[j] = -(1 + 2 * self.sigma)
                self.c[j] = self.sigma
                self.d[j] = -2 * self.u[k - 1][j] + self.u[k - 2][j]

            if self.type == '1-2':
                self.b[0] = self.alpha / self.h / (self.beta - self.alpha / self.h)
                self.c[0] = 1
                self.d[0] = 1 / (self.beta - self.alpha / self.h) * self.phi0(k * self.tau)
                self.a[-1] = -self.gamma / self.h / (self.delta + self.gamma / self.h)
                self.d[-1] = 1 / (self.delta + self.gamma / self.h) * self.phi_l(k * self.tau)

            elif self.type == '2-2':
                self.b[0] = 2 * self.a_const / self.h
                self.c[0] = -2 * self.a_const / self.h + self.h / self.tau ** 2 - self.c_const * self.h + -self.d_const * self.h / (2 * self.tau) + self.beta / self.alpha * (2 * self.a_const + self.b_const * self.h)
                self.d[0] = self.h / self.tau ** 2 * (self.u[k - 2][0] - 2 * self.u[k - 1][0]) - self.h * self.function() + -self.d_const * self.h / (2 * self.tau) * self.u[k - 2][0] + (2 * self.a_const - self.b_const * self.h) / self.alpha * self.phi0(k * self.tau)
                self.a[-1] =-self.b[0]
                self.d[-1] = self.h / self.tau ** 2 * (-self.u[k - 2][0] + 2 * self.u[k - 1][0]) + self.h * self.function() + self.d_const * self.h / (2 * self.tau) * self.u[k - 2][0] + (2 * self.a_const + self.b_const * self.h) / self.alpha * self.phi_l(k * self.tau)

            elif self.type == '2-3':
                k1 = 2 * self.h * self.beta - 3 * self.alpha
                k2 = 2 * self.h * self.delt + 3 * self.gamma
                omega = self.tau ** 2 * self.b_const / (2 * self.h)
                xi = self.d_const * self.tau / 2
                self.b[0] = 4 * self.alpha - self.alpha / (self.sigma + omega) * (1 + xi + 2 * self.sigma - self.c_const * self.tau ** 2)
                self.c[0] = k1 - self.alpha * (omega - self.sigma) / (omega + self.sigma)
                self.d[0] = 2 * self.h * self.phi0(k * self.tau) + self.alpha * self.d[1] / (-self.sigma - omega)
                self.a[-1] = -self.gamma / (omega - self.sigma) * (1 + xi + 2 * self.sigma - self.c_const * self.tau ** 2) - 4 * self.gamma
                self.d[-1] = 2 * self.h * self.phi_l(k * self.tau) - self.gamma * self.d[-2] / (omega - self.sigma)

            self.u[k] = self.run_through_method()

        return self.u

    def Explicit(self, N, K, T):
        self.u = np.zeros((K, N))
        for j in range(N):
            x_current = j * self.h
            self.u[0][j] = self.psi1(x_current)
            if self.accuracy_level == 1:
                self.u[1][j] = self.psi1(x_current) + self.psi2(x_current) * self.tau + self.psi12(x_current) * self.tau ** 2 / 2
            elif self.accuracy_level == 2:
                k = self.tau ** 2 / 2
                self.u[1][j] = (1 + self.c_const * k) * self.psi2(x_current) + self.a_const * k * self.psi12(x_current) + self.b_const * k * self.psi11(x_current) + (self.tau - self.d_const * k) * self.psi1(x_current) + k * self.function()

        if self.type == '1-2':
            LBound = self.l_bound_12
            RBound = self.r_bound_12
        elif self.type == '2-2':
            LBound = self.l_bound_22
            RBound = self.r_bound_22
        elif self.type == '2-3':
            LBound = self.l_bound_23
            RBound = self.r_bound_23
        for k in range(2, K):
            t = k * self.tau
            for j in range(1, N - 1):
                self.u[k][j] = self.u[k - 1][j + 1] * (self.sigma + self.b_const * self.tau ** 2 / (2 * self.h)) + self.u[k - 1][j] * (-2 * self.sigma + 2 + self.c_const * self.tau ** 2) + self.u[k - 1][j - 1] * (self.sigma - self.b_const * self.tau ** 2 / (2 * self.h)) - self.u[k - 2][j] + self.tau ** 2 * self.function()

            self.u[k][0] = LBound(k, t)
            self.u[k][-1] = RBound(k, t)

        return self.u

    def l_bound_12(self, k, t):
        return -(self.alpha / self.h) / (self.beta - self.alpha / self.h) * self.u[k - 1][1] + self.phi0(t) / (self.beta - self.alpha / self.h)

    def r_bound_12(self, k, t):
        return (self.gamma / self.h) / (self.delta + self.gamma / self.h) * self.u[k - 1][-2] + self.phi_l(t) / (self.delta + self.gamma / self.h)

    def l_bound_22(self, k, t):
        n = self.c_const * self.h - 2 * self.a_const / self.h - self.h / self.tau ** 2 - self.d_const * self.h / (2 * self.tau) + self.beta / self.alpha * (2 * self.a_const - self.b_const * self.h)
        return 1 / n * (- 2 * self.a_const / self.h * self.u[k][1] + self.h / self.tau ** 2 * (self.u[k - 2][0] - 2 * self.u[k - 1][0]) + -self.d_const * self.h / (2 * self.tau) * self.u[k - 2][0] + -self.h * self.function() + (2 * self.a_const - self.b_const * self.h) / self.alpha * self.phi0(t))

    def r_bound_22(self, k, t):
        n = -self.c_const * self.h + 2 * self.a_const / self.h + self.h / self.tau ** 2 + self.d_const * self.h / (2 * self.tau) + self.delta / self.gamma * (2 * self.a_const + self.b_const * self.h)
        return 1 / n * (2 * self.a_const / self.h * self.u[k][-2] + self.h / self.tau ** 2 * (2 * self.u[k - 1][-1] - self.u[k - 2][-1]) + self.d_const * self.h / (2 * self.tau) * self.u[k - 2][-1] + self.h * self.function() + (2 * self.a_const + self.b_const * self.h) / self.gamma * self.phi_l(t))

    def l_bound_23(self, k, t):
        n = 2 * self.h * self.beta - 3 * self.alpha
        return self.alpha / n * self.u[k - 1][2] - 4 * self.alpha / n * self.u[k - 1][1] + 2 * self.h / n * self.phi0(t)

    def r_bound_23(self, k, t):
        n = 2 * self.h * self.delta + 3 * self.gamma
        return 4 * self.gamma / n * self.u[k - 1][-2] - self.gamma / n * self.u[k - 1][-3] + 2 * self.h / n * self.phi_l(t)