import matplotlib.pyplot as plt
import copy
import numpy as np


def calculate_difference(L, u, nx, ny):
    mx = 0
    for i in range(nx):
        for j in range(ny):
            mx = max(mx, abs(u[i][j] - L[i][j]))
    return mx


class ProblemData:
    def __init__(self, parameters):
        self.a = parameters['a']
        self.b = parameters['b']
        self.c = parameters['c']
        self.d = parameters['d']
        self.lx = parameters['lx']
        self.ly = parameters['ly']
        self.w = parameters['w']
        self.f = parameters['f']
        self.alpha1 = parameters['alpha1']
        self.alpha2 = parameters['alpha2']
        self.beta1 = parameters['beta1']
        self.beta2 = parameters['beta2']
        self.gamma1 = parameters['gamma1']
        self.gamma2 = parameters['gamma2']
        self.delta1 = parameters['delta1']
        self.delta2 = parameters['delta2']
        self.phi1 = parameters['phi1']
        self.phi2 = parameters['phi2']
        self.phi3 = parameters['phi3']
        self.phi4 = parameters['phi4']
        self.solution = parameters['solution']


class EllipticalSolver:
    def __init__(self, parameters, equation_type):
        self.data = ProblemData(parameters)
        self.hx = 0
        self.hy = 0
        self.iteration = 0
        self.eps = 1e-6
        try:
            self.solve_function = getattr(self, f'{equation_type}_solver')
        except:
            raise Exception("This type does not exist")

    def initialize_u(self, x, y):
        u = np.zeros((len(x), len(y)))
        for i in range(len(x)):
            u[i][0] = self.data.phi3(x[i]) / self.data.gamma2
            u[i][-1] = self.data.phi4(x[i]) / self.data.delta2
        for j in range(len(y)):
            u[0][j] = self.data.phi1(y[j]) / self.data.alpha2
            u[-1][j] = self.data.phi2(y[j]) / self.data.beta2

        return u

    def solve(self, nx, ny):
        self.hx = self.data.lx / nx
        self.hy = self.data.ly / ny
        x = np.arange(0, self.data.lx + self.hx, self.hx)
        y = np.arange(0, self.data.ly + self.hy, self.hy)

        u = self.initialize_u(x, y)
        for i in range(1, nx):
            for j in range(1, ny):
                u[i][j] = u[0][j] + (x[i] - x[0]) * (u[-1][j] - u[0][j]) / (x[-1] - x[0])

        return self.solve_function(nx, ny, x, y, u)

    def analytic_solve(self, nx, ny):
        self.hx = self.data.lx / nx
        self.hy = self.data.ly / ny
        x = np.arange(0, self.data.lx + self.hx, self.hx)
        y = np.arange(0, self.data.ly + self.hy, self.hy)

        u = []
        for yi in y:
            u.append([self.data.solution(xi, yi) for xi in x])
        return u

    def simple_iteration_method_solver(self, nx, ny, x, y, u):
        cur_eps = 1e9
        while self.iteration < 10000:
            L = copy.deepcopy(u)
            u = self.initialize_u(x, y)
            for j in range(1, len(y) - 1):
                for i in range(1, len(x) - 1):
                    u[i][j] = (self.hx * self.hx * self.data.f(x[i], y[j]) -
                               (L[i + 1][j] + L[i - 1][j]) - self.data.d * self.hx * self.hx *
                               (L[i][j + 1] + L[i][j - 1]) /
                               (self.hy * self.hy) - self.data.a * self.hx * 0.5 *
                               (L[i + 1][j] - L[i - 1][j]) - self.data.b * self.hx * self.hx *
                               (L[i][j + 1] - L[i][j - 1]) /
                               (2 * self.hy)) / (self.data.c * self.hx * self.hx - 2 *
                                                 (self.hy * self.hy + self.data.d * self.hx * self.hx) /
                                                 (self.hy * self.hy))
            last_eps = cur_eps
            cur_eps = calculate_difference(L, u, nx, ny)
            if calculate_difference(L, u, nx, ny) <= self.eps or last_eps < cur_eps:
                break
            self.iteration += 1
        return u, self.iteration

    def seidel_method_solver(self, nx, ny, x, y, u):
        cur_eps = 1e9
        while self.iteration < 10000:
            L = copy.deepcopy(u)
            u = self.initialize_u(x, y)
            for j in range(1, len(y) - 1):
                for i in range(1, len(x) - 1):
                    u[i][j] = ((self.hx ** 2) * self.data.f(x[i], y[j]) -
                               (L[i + 1][j] + u[i - 1][j]) - self.data.d * (self.hx ** 2) *
                               (L[i][j + 1] + u[i][j - 1]) / (self.hy ** 2) - self.data.a * self.hx * 0.5 *
                               (L[i + 1][j] - u[i - 1][j]) - self.data.b * (self.hx ** 2) *
                               (L[i][j + 1] - u[i][j - 1]) /
                               (2 * self.hy)) / \
                              (self.data.c * (self.hx ** 2) - 2 * (self.hy ** 2 + self.data.d * (self.hx ** 2)) /
                               (self.hy ** 2))
            last_eps = cur_eps
            cur_eps = calculate_difference(L, u, nx, ny)
            if cur_eps <= self.eps or last_eps < cur_eps:
                break
            self.iteration += 1
        return u, self.iteration

    def simple_iteration_method_relaxed_solver(self, nx, ny, x, y, u):
        cur_eps = 1e9
        while self.iteration < 10000:
            L = copy.deepcopy(u)
            u = self.initialize_u(x, y)
            for j in range(1, len(y) - 1):
                for i in range(1, len(x) - 1):
                    u[i][j] = (((self.hx ** 2) * self.data.f(x[i], y[j]) -
                                (L[i + 1][j] + u[i - 1][j]) - self.data.d * (self.hx ** 2) *
                                (L[i][j + 1] + u[i][j - 1]) / (self.hy ** 2) - self.data.a * self.hx * 0.5 *
                                (L[i + 1][j] - u[i - 1][j]) - self.data.b * (self.hx ** 2) *
                                (L[i][j + 1] - u[i][j - 1]) /
                                (2 * self.hy)) / (self.data.c * (self.hx ** 2) - 2 *
                                                  (self.hy ** 2 + self.data.d * (self.hx ** 2)) /
                                                  (self.hy ** 2))) * self.data.w + (1 - self.data.w) * L[i][j]
            last_eps = cur_eps
            cur_eps = calculate_difference(L, u, nx, ny)
            if calculate_difference(L, u, nx, ny) <= self.eps or last_eps < cur_eps:
                break
            self.iteration += 1
        return u, self.iteration


def calculate_error(a, b):
    err = 0
    lst = [abs(i - j) for i, j in zip(a, b)]
    for each in lst:
        err = max(err, each)
    return err


parameters = {'equation_type': 'simple_iteration_method_relaxed', 'nx': 40, 'ny': 40}
