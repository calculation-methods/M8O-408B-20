import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

plt.rcParams['figure.figsize'] = [8, 7]

mu1 = 1
mu2 = 1
a = 1

Nx = 50
Ny = 50
Nt = 40
lx = mu1 * (np.pi / 2)
ly = mu2 * (np.pi / 2)
T = 2
hx = lx / Nx
hy = ly / Ny
tau = T / Nt

def phi1(x, t):
    return np.cos(mu1 * x) * np.exp(-(mu1 * mu1 + mu2 * mu2) * a * t)

def phi2(x, t):
    return 0

def phi3(y, t):
    return np.cos(mu2 * y) * np.exp(-(mu1 * mu1 + mu2 * mu2) * a * t)

def phi4(y, t):
    return 0

def psi(x, y):
    return np.cos(mu1 * x) * np.cos(mu2 * y)

def U(x, y, t):
    return np.cos(mu1 * x) * np.cos(mu2 * y) * np.exp(-(mu1 * mu1 + mu2 * mu2) * a * t)

def norm(v1, v2):
    return np.amax(np.abs(v1 - v2))

def Check(A):
    if np.shape(A)[0] != np.shape(A)[1]:
        return False
    n = np.shape(A)[0]
    for i in range(n):
        sum = 0
        for j in range(n):
            if i != j:
                sum += abs(A[i][j])
        if abs(A[i][i]) < sum:
            return False
    return True

def solve(a, b):
    if (Check(a)):
        p = np.zeros(len(b))
        q = np.zeros(len(b))
        p[0] = -a[0][1] / a[0][0]
        q[0] = b[0] / a[0][0]
        for i in range(1, len(p) - 1):
            p[i] = -a[i][i + 1] / (a[i][i] + a[i][i - 1] * p[i - 1])
            q[i] = (b[i] - a[i][i - 1] * q[i - 1]) / (a[i][i] + a[i][i - 1] * p[i - 1])
        i = len(a) - 1
        p[-1] = 0
        q[-1] = (b[-1] - a[-1][-2] * q[-2]) / (a[-1][-1] + a[-1][-2] * p[-2])
        x = np.zeros(len(b))
        x[-1] = q[-1]
        for i in reversed(range(len(b) - 1)):
            x[i] = p[i] * x[i + 1] + q[i]
        return x

def error(Nt, l, tau, U):
    NArray = [10, 20, 40]
    size = np.size(NArray)
    hArray = np.zeros(size)
    errors1x = np.zeros(size)
    errors2x = np.zeros(size)
    errors1y = np.zeros(size)
    errors2y = np.zeros(size)
    for i in range(0, size):
        hArray[i] = l / NArray[i]
        u1 = VariableDirectionMethod(Nt, NArray[i], NArray[i], tau, hArray[i], hArray[i])
        u2 = FractionalStepsMethod(Nt, NArray[i], NArray[i], tau, hArray[i], hArray[i])
        t = tau * Nt / 2
        x = hArray[i] * NArray[i] / 2
        y = hArray[i] * NArray[i] / 2
        uxCorrect = np.array([U(xi * hArray[i], y, t) for xi in range(NArray[i] + 1)])
        uyCorrect = np.array([U(x, yi * hArray[i], t) for yi in range(NArray[i] + 1)])
        u1xCalc = u1[int(Nt / 2)][:][int(NArray[i] / 2)]
        u2xCalc = u2[int(Nt / 2)][:][int(NArray[i] / 2)]
        u1yCalc = u1[int(Nt / 2)][int(NArray[i] / 2)][:]
        u2yCalc = u2[int(Nt / 2)][int(NArray[i] / 2)][:]
        errors1x[i] = np.amax(np.abs(uxCorrect - u1xCalc))
        errors2x[i] = np.amax(np.abs(uxCorrect - u2xCalc))
        errors1y[i] = np.amax(np.abs(uyCorrect - u1yCalc))
        errors2y[i] = np.amax(np.abs(uyCorrect - u2yCalc))
    return NArray, errors1x, errors2x, errors1y, errors2y

def showErrors(Nt, l, tau, U):
    NArray, errors1x, errors2x, errors1y, errors2y = error(Nt, l, tau, U)
    colors = ['blue', 'red']
    delta = np.zeros(np.size(NArray))
    for i in range(np.size(NArray)):
        delta[i] = l / NArray[i]
    delta2 = np.zeros(np.size(NArray))
    for i in range(np.size(NArray)):
        delta2[i] = l / NArray[np.size(NArray) - i - 1]
    fig, ax = plt.subplots()
    plt.plot(delta, errors1x, color=colors[0], label='Метод переменныхнаправлений')
    plt.plot(delta2, errors2x, color=colors[1], label='Метод дробных шагов')
    ax.set_xlabel('delta X')
    ax.set_ylabel('Epsilon')
    plt.grid()
    ax.legend()
    plt.show()
    fig, ax = plt.subplots()
    plt.plot(delta, errors1y, color=colors[0], label='Метод переменныхнаправлений')
    plt.plot(delta2, errors2y, color=colors[1], label='Метод дробных шагов')
    ax.set_xlabel('delta Y')
    ax.set_ylabel('Epsilon')
    plt.grid()
    ax.legend()
    plt.show()

def showSolution(Nx, Ny, Nt, hx, hy, tau, U, u):
    xArray = np.array([i * hx for i in range(Nx + 1)])
    yArray = np.array([j * hy for j in range(Ny + 1)])
    fig, ax = plt.subplots(2)
    t = [int(Nt * 0.05), int(Nt * 0.4), int(Nt * 0.7)]
    xFix = int(Nx / 2)
    yFix = int(Ny / 4)
    colors = ['blue', 'green', 'red']
    for i in range(len(t)):
        uCorrect = np.zeros(Nx + 1)
        for x in range(Nx + 1):
            uCorrect[x] = U(x * hx, yFix * hy, t[i] * tau)
        uCalc = u[t[i]][:][yFix]
        ax[0].plot(yArray, uCorrect, color=colors[i], label='t=%s' % round(t[i] * tau, 2))
        ax[0].plot(yArray, uCalc, color=colors[i], linestyle='--')
    for i in range(len(t)):
        uCorrect = np.zeros(Ny + 1)
        for y in range(Ny + 1):
            uCorrect[y] = U(xFix * hx, y * hy, t[i] * tau)
        uCalc = u[t[i]][xFix][:]
        ax[1].plot(yArray, uCorrect, color=colors[i], label='t=%s' % round(t[i] * tau, 2))
        ax[1].plot(yArray, uCalc, color=colors[i], linestyle='--')
    label1 = 'x (yFix=%s)' % round(yFix * hy, 2)
    label2 = 'y (xFix=%s)' % round(xFix * hx, 2)
    ax[0].set_xlabel(label1)
    ax[0].set_ylabel('U(x, y, t)')
    ax[0].grid()
    ax[0].legend()
    ax[1].set_xlabel(label2)
    ax[1].set_ylabel('U(x, y, t)')
    ax[1].grid()
    ax[1].legend()
    plt.show()

def VariableDirectionMethod(Nt, Nx, Ny, tau, hx, hy):
    u = np.zeros((Nt + 1, Nx + 1, Ny + 1))
    for t in range(Nt + 1):
        for x in range(Nx + 1):
            u[t][x][0] = phi1(x * hx, t * tau)
            u[t][x][Ny] = phi2(x * hx, t * tau)
    for t in range(Nt + 1):
        for y in range(Ny + 1):
            u[t][0][y] = phi3(y * hy, t * tau)
            u[t][Nx][y] = phi4(y * hy, t * tau)
    for x in range(Nx + 1):
        for y in range(Ny + 1):
            u[0][x][y] = psi(y * hy, x * hx)
    for t in range(Nt):
        tmp = deepcopy(u[t])
    for y in range(1, Ny):
        matrix = np.zeros((Nx - 1, Nx - 1))
        d = np.zeros(Nx - 1)
        ai = a * tau / (2 * hx * hx)
        bi = -(a * tau / (hx * hx) + 1)
        ci = a * tau / (2 * hx * hx)
        matrix[0][0] = bi
        matrix[0][1] = ci
        d[0] = -(u[t][1][y] + (a * tau / (2 * hy * hy)) * (u[t][1][y - 1] - 2 * u[t][1][y] + u[t][1][y + 1]) + a * tau / (2 * hx * hx) * phi3(y * hy, (t + 1 / 2) * tau))
        for x in range(1, Nx - 2):
            matrix[x][x - 1] = ai
            matrix[x][x] = bi
            matrix[x][x + 1] = ci
            d[x] = -(u[t][x + 1][y] + (a * tau / (2 * hy * hy)) * (u[t][x + 1][y - 1] - 2 * u[t][x + 1][y] + u[t][x + 1][y + 1]))
        matrix[Nx - 2][Nx - 3] = ai
        matrix[Nx - 2][Nx - 2] = bi
        d[Nx - 2] = -(u[t][Nx - 1][y] + (a * tau / (2 * hy * hy)) * (u[t][Nx - 1][y - 1] - 2 * u[t][Nx - 1][y] + u[t][Nx - 1][y + 1]) + a * tau / (2 * hx * hx) * phi4(y * hy, (t + 1 / 2) * tau))
        ans = np.linalg.solve(matrix, d)
        tmp[1:Nx, y] = ans
    tmp[0][:] = np.array([phi3(j * hy, (t + 1 / 2) * tau) for j in range(Ny + 1)])
    tmp[Nx][:] = np.array([phi4(j * hy, (t + 1 / 2) * tau) for j in range(Ny + 1)])
    tmp[:][0] = np.array([phi1(i * hx, (t + 1 / 2) * tau) for i in range(Nx + 1)])
    tmp[:][Ny] = np.array([phi2(i * hx, (t + 1 / 2) * tau) for i in range(Nx + 1)])
    for x in range(1, Nx):
        matrix = np.zeros((Ny - 1, Ny - 1))
        d = np.zeros(Ny - 1)
        ai = a * tau / (2 * hy * hy)
        bi = -(a * tau / (hy * hy) + 1)
        ci = a * tau / (2 * hy * hy)
        matrix[0][0] = bi
        matrix[0][1] = ci
        d[0] = -(tmp[x][1] + (a * tau / (2 * hx * hx)) * (tmp[x - 1][1] - 2 * tmp[x][1] + tmp[x + 1][1]) + a * tau / (2 * hy * hy) * phi1(x * hx, (t + 1) * tau))
        for y in range(1, Ny - 2):
            matrix[y][y - 1] = ai
            matrix[y][y] = bi
            matrix[y][y + 1] = ci
            d[y] = -(tmp[x][y + 1] + (a * tau / (2 * hx * hx)) * (tmp[x - 1][y + 1] - 2 * tmp[x][y + 1] + tmp[x + 1][y + 1]))
        matrix[Ny - 2][Ny - 3] = ai
        matrix[Ny - 2][Ny - 2] = bi
        d[Ny - 2] = -(tmp[x][Ny - 1] + (a * tau / (2 * hx * hx)) * (tmp[x - 1][Ny - 1] - 2 * tmp[x][Ny - 1] + tmp[x + 1][Ny - 1]) + a * tau / (2 * hy * hy) * phi2(x * hx, (t + 1) * tau))
        ans = np.linalg.solve(matrix, d)
        u[t + 1, x, 1:Ny] = ans
    return u


def FractionalStepsMethod(Nt, Nx, Ny, tau, hx, hy):
    u = np.zeros((Nt + 1, Nx + 1, Ny + 1))
    for t in range(Nt + 1):
        for x in range(Nx + 1):
            u[t][x][0] = phi1(x * hx, t * tau)
            u[t][x][Ny] = phi2(x * hx, t * tau)
    for t in range(Nt + 1):
        for y in range(Ny + 1):
            u[t][0][y] = phi3(y * hy, t * tau)
            u[t][Nx][y] = phi4(y * hy, t * tau)
    for x in range(Nx + 1):
        for y in range(Ny + 1):
            u[0][x][y] = psi(y * hy, x * hx)
    for t in range(Nt):
        tmp = deepcopy(u[t])
        for y in range(1, Ny):
            matrix = np.zeros((Nx - 1, Nx - 1))
            d = np.zeros(Nx - 1)
            ai = a * tau / (hx * hx)
            bi = -(a * tau * 2 / (hx * hx) + 1)
            ci = a * tau / (hx * hx)
            matrix[0][0] = bi
            matrix[0][1] = ci
            d[0] = -(u[t][1][y] + a * tau / (hx * hx) * phi3(y * hy, (t + 1 / 2) * tau))
            for x in range(1, Nx - 2):
                matrix[x][x - 1] = ai
                matrix[x][x] = bi
                matrix[x][x + 1] = ci
                d[x] = -u[t][x + 1][y]
            matrix[Nx - 2][Nx - 3] = ai
            matrix[Nx - 2][Nx - 2] = bi
            d[Nx - 2] = -(u[t][Nx - 1][y] + a * tau / (hx * hx) * phi4(y * hy, (t + 1 / 2) * tau))
            ans = np.linalg.solve(matrix, d)
            tmp[1:Nx, y] = ans
        tmp[0][:] = np.array([phi3(j * hy, (t + 1 / 2) * tau) for j in range(Ny + 1)])
        tmp[Nx][:] = np.array([phi4(j * hy, (t + 1 / 2) * tau) for j in range(Ny + 1)])
        tmp[:][0] = np.array([phi1(i * hx, (t + 1 / 2) * tau) for i in range(Nx + 1)])
        tmp[:][Ny] = np.array([phi2(i * hx, (t + 1 / 2) * tau) for i in range(Nx + 1)])
        for x in range(1, Nx):
            matrix = np.zeros((Ny - 1, Ny - 1))
            d = np.zeros(Ny - 1)
            ai = a * tau / (hy * hy)
            bi = -(a * tau * 2 / (hy * hy) + 1)
            ci = a * tau / (hy * hy)
            matrix[0][0] = bi
            matrix[0][1] = ci
            d[0] = -(tmp[x][1] + a * tau / (hy * hy) * phi1(x * hx, (t + 1) * tau))
            for y in range(1, Ny - 2):
                matrix[y][y - 1] = ai
                matrix[y][y] = bi
                matrix[y][y + 1] = ci
                d[y] = -tmp[x][y + 1]
            matrix[Ny - 2][Ny - 3] = ai
            matrix[Ny - 2][Ny - 2] = bi
            d[Ny - 2] = -(tmp[x][Ny - 1] + a * tau / (hy * hy) * phi2(x * hx, (t + 1) * tau))
            ans = solve(matrix, d)
            u[t + 1, x, 1:Ny] = ans
    return u

def main():
    u1 = VariableDirectionMethod(Nt, Nx, Ny, tau, hx, hy)
    showSolution(Nx, Ny, Nt, hx, hy, tau, U, u1)
    u2 = FractionalStepsMethod(Nt, Nx, Ny, tau, hx, hy)
    showSolution(Nx, Ny, Nt, hx, hy, tau, U, u2)
    showErrors(Nt, lx, tau, U)


main()
