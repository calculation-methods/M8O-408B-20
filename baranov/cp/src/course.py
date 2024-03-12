import numpy as np
from numpy.linalg import norm
from scipy.sparse import diags, csc_matrix
from time import time

def getMatrix(filename):
    with open(filename) as f:
        shape = int(f.readline())
        matrix = [[float(num) for num in line.split()]
                  for _, line in zip(range(shape), f)]
        matrix = csc_matrix(matrix)
        b = np.array([float(num) for num in f.readline().split()])
        return matrix, b

def biCGStabSolve(matrix, b, eps, shape, x0, k):
    r0 = b - matrix @ x0
    x0 = x0
    r2 = r0
    rho0 = 1
    alpha0 = 1
    omega0 = 1
    v0 = np.array([0] * shape)
    p0 = np.array([0] * shape)
    while True:
        rho = r2 @ r0
        beta = (rho * alpha0) / (rho0 * omega0)
        p = r0 + beta * (p0 - omega0 * v0)
        v = matrix @ p
        alpha = rho / (r2 @ v)
        s = r0 - alpha * v
        t = matrix @ s
        omega = (t @ s) / (t @ t)
        x = x0 + omega * s + alpha * p
        r = s - omega * t

        k += 1
        if norm(r) < eps:
            break
        r0 = r
        rho0 = rho
        alpha0 = alpha
        omega0 = omega
        v0 = v
        p0 = p
        x0 = x
    return x


def printSolution(matrix, b):
    eps = 1e-5
    shape = matrix.shape[0]
    x0 = np.array([0] * shape)
    k = 0

    start = time()
    x = biCGStabSolve(matrix, b, eps, shape, x0, k)
    end = time()
    start2 = time()
    x2 = np.linalg.solve(matrix.toarray(), b)
    end2 = time()
    print('My solve:\n')
    print(f'{x.round(5)}\n')
    print(f'EPS = {eps}\n')
    print(f'Shape = {shape}\n')
    print(f'Count of iterations = {k}\n')
    print(f'Mean = {np.mean(x)}\n')
    print(f'Time = {round(end - start, 5)} sec\n')
    print('\nNumPy solve:\n')
    print(f'{x2.round(5)}\n')
    print(f'Mean = {np.mean(x2)}\n')
    print(f'Time = {round(end2 - start2, 5)} sec\n')


for fileName in ['test10', 'test20', 'test30']:
    matrix, b = getMatrix('test10')
    printSolution(matrix, b)