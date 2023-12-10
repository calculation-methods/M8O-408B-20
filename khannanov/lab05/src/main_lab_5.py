import csv
import math
import numpy as np


class ParabolicSolver:
    def __init__(self, L, psi, function, phi0, phiL, analytical_solution, ptype):
        self.L = L
        self.psi = psi
        self.function = function
        self.phi0 = phi0
        self.phiL = phiL
        self.analytical_solution = analytical_solution
        self.type = ptype

    def implicit(self, N, K, T):
        h = self.L / N
        tau = T / K
        sigma = tau / (h * h)
        u = np.zeros((K, N))

        for i in range(1, N - 1):
            u[0, i] = self.psi(i * h)
        u[0, -1] = 0

        for k in range(1, K):
            a = np.zeros(N)
            b = np.zeros(N)
            c = np.zeros(N)
            d = np.zeros(N)

            for j in range(1, N - 1):
                a[j] = sigma
                b[j] = -(1 + 2 * sigma)
                c[j] = sigma
                d[j] = -u[k - 1, j] - tau * self.function(j * h, k * tau)

            a[0] = 0
            b[0] = -(1 + 2 * sigma)
            c[0] = sigma
            a[-1] = sigma
            b[-1] = -(1 + 2 * sigma)
            c[-1] = 0

            if self.type == "1-1":
                d[0] = -(u[k - 1, 0] + sigma * self.phi0(k * tau))
                d[-1] = -(u[k - 1, -1] + sigma * self.phiL(k * tau))
            elif self.type == "2-2":
                d[0] = -(u[k - 1, 0] + sigma * self.phi0(k * tau) - tau * self.function(0, k * tau))
                d[-1] = -(u[k - 1, -1] + sigma * self.phiL(k * tau) - tau * self.function((N - 1) * h, k * tau))
            elif self.type == "2-3":
                d[0] = -((1 - sigma) * u[k - 1, 1] + sigma / 2 * u[k - 1, 0] - tau * self.function(0, k * tau) - sigma * self.phi0(k * tau))
                d[-1] = self.phiL(k * tau)
                d[-1] += (tau * self.function((N - 1) * h, k * tau) * h) / (2 * tau * u[k - 1, -1] + 1)

            for i in range(1, N):
                m = a[i] / b[i - 1]
                b[i] -= m * c[i - 1]
                d[i] -= m * d[i - 1]

            u[k, -1] = d[-1] / b[-1]
            for i in range(N - 2, -1, -1):
                u[k, i] = (d[i] - c[i] * u[k, i + 1]) / b[i]

        return u

    def explicit(self, N, K, T):
        h = self.L / N
        tau = T / K
        sigma = tau / (h * h)
        u = np.zeros((K, N))

        for j in range(1, N - 1):
            u[0, j] = self.psi(j * h)

        for k in range(1, K):
            u[k, 0] = self.phi0(k * tau)
            for j in range(1, N - 1):
                u[k, j] = sigma * u[k - 1, j + 1] + (1 - 2 * sigma) * u[k - 1, j] + sigma * u[k - 1, j - 1] + tau * self.function(j * h, k * tau)

            if self.type == "1-1":
                u[k, -1] = u[k, -2] + self.phiL(k * tau) * h
            elif self.type == "2-2":
                u[k, -1] = self.phiL(k * tau)
            elif self.type == "2-3":
                u[k, -1] = (self.phiL(k * tau) + u[k, -2] / h + 2 * tau * u[k - 1, -1] / h) / (1 / h + 2 * tau / h)

        return u

    def crank_nicholson(self, N, K, T):
        h = self.L / N
        tau = T / K
        sigma = tau / (h * h)
        u = np.zeros((K, N))

        for j in range(1, N - 1):
            u[0, j] = self.psi(j * h)

        for k in range(1, K):
            a = np.zeros(N)
            b = np.zeros(N)
            c = np.zeros(N)
            d = np.zeros(N)
            u_implicit = np.zeros(N)

            for j in range(1, N - 1):
                a[j] = sigma
                b[j] = -(1 + 2 * sigma)
                c[j] = sigma
                d[j] = -u[k - 1, j] - tau * self.function(j * h, k * tau)

            a[0] = 0
            b[0] = -(1 + 2 * sigma)
            c[0] = sigma
            a[-1] = sigma
            b[-1] = -(1 + 2 * sigma)
            c[-1] = 0

            if self.type == "1-1":
                d[0] = -(u[k - 1, 0] + sigma * self.phi0(k * tau))
                d[-1] = -(u[k - 1, -1] + sigma * self.phiL(k * tau))
            elif self.type == "2-2":
                d[0] = -(u[k - 1, 0] + sigma * self.phi0(k * tau) - tau * self.function(0, k * tau))
                d[-1] = -(u[k - 1, -1] + sigma * self.phiL(k * tau) - tau * self.function((N - 1) * h, k * tau))
            elif self.type == "2-3":
                d[0] = -((1 - sigma) * u[k - 1, 1] + sigma / 2 * u[k - 1, 0] - tau * self.function(0, k * tau) - sigma * self.phi0(k * tau))
                d[-1] = self.phiL(k * tau)
                d[-1] += (tau * self.function((N - 1) * h, k * tau) * h) / (2 * tau * u[k - 1, -1] + 1)

            p = np.zeros(N)
            q = np.zeros(N)
            p[0] = -c[0] / b[0]
            q[0] = d[0] / b[0]

            for i in range(1, N):
                m = a[i] / b[i - 1]
                b[i] -= m * c[i - 1]
                d[i] -= m * d[i - 1]

            u_implicit[-1] = q[-1]
            for i in range(N - 2, -1, -1):
                u_implicit[i] = p[i] * u_implicit[i + 1] + q[i]

            u_explicit = np.zeros(N)
            for j in range(1, N - 1):
                u_explicit[j] = sigma * u[k - 1, j + 1] + (1 - 2 * sigma) * u[k - 1, j] + sigma * u[k - 1, j - 1] + tau * self.function(j * h, k * tau)

            if self.type == "1-1":
                u_implicit[-1] = u[k, -2] + self.phiL(k * tau) * h
            elif self.type == "2-2":
                u_implicit[-1] = self.phiL(k * tau)
            elif self.type == "2-3":
                u_implicit[-1] = (self.phiL(k * tau) + u[k, -2] / h + 2 * tau * u[k - 1, -1] / h) / (1 / h + 2 * tau / h)

            for j in range(N):
                u[k, j] = 0.5 * u_implicit[j] + 0.5 * u_explicit[j]

        return u


def analytical_solution_matrix(N, K, T, solver):
    h = solver.L / N
    tau = T / K
    u = np.zeros((K, N))

    for k in range(K):
        for j in range(N):
            u[k, j] = solver.analytical_solution(j * h, k * tau)

    return u
