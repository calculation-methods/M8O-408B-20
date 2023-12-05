import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from objects import HyperbolicSolver

def plot(dict_, N, K, T, save_file = "borisov/lab6/plot.png"):
    fig = plt.figure(figsize = plt.figaspect(0.7))

    x = np.arange(0, np.pi, np.pi / N)
    t = np.arange(0, T, T / K)
    x, t = np.meshgrid(x, t)
    z1 = np.array(dict_['numerical'])
    z2 = np.array(dict_['analytic'])

    ax = fig.add_subplot(1, 2, 1, projection = '3d')
    ax.set_xlabel('x', fontsize = 20)
    ax.set_ylabel('t', fontsize = 20)
    ax.set_zlabel('u', fontsize = 20)
    plt.title('numerical')
    surf = ax.plot_surface(x, t, z1, cmap = cm.coolwarm, linewidth = 0, antialiased = True)
    fig.colorbar(surf, shrink = 0.5, aspect = 15)

    ax = fig.add_subplot(1, 2, 2, projection = '3d')
    ax.set_xlabel('x', fontsize = 20)
    ax.set_ylabel('t', fontsize = 20)
    ax.set_zlabel('u', fontsize = 20)
    plt.title('analytic')
    surf = ax.plot_surface(x, t, z2, cmap = cm.coolwarm, linewidth = 0, antialiased = True)

    fig.colorbar(surf, shrink = 0.5, aspect = 15)

    plt.savefig(save_file)
    plt.show()


def compare_error(dict_):
    error = [[abs(i - j) for i, j in zip(x, y)] for x, y in zip(dict_['numerical'], dict_['analytic'])]
    return error

data = {'equation_type': 'explicit', 'N': 50, 'K': 100, 'T': 1}

if __name__ == '__main__':
    equation_type = data['equation_type']
    N, K, T = int(data['N']), int(data['K']), int(data['T'])

    params = {
        'a': 1,
        'b': 0,
        'c': -3,
        'd': 0,
        'l': np.pi,
        'f': lambda: 0,
        'alpha': 1,
        'beta': 0,
        'gamma': 1,
        'delta': 0,
        'psi1': lambda x: 0,
        'psi2': lambda x: 2 * np.cos(x),
        'psi1_dir1': lambda x: 0,
        'psi1_dir2': lambda x: -2 * np.sin(x),
        'phi0': lambda t: np.sin(2 * t),
        'phi1': lambda t: -1 * np.sin(2 * t),
        'bound_type': 'a1p2',
        'approximation': 'p1',
        'solution': lambda x, t: np.cos(x) * np.sin(2 * t),
    }

    var = HyperbolicSolver(params, equation_type)

    ans = {
        'numerical': var.solve(N, K, T).tolist(),
        'analytic': var.analyticSolve(N, K, T).tolist()
    }

    plot(ans, N, K, T)

    error = compare_error(ans)
    avg_err = 0.0
    for i in error:
        for j in i:
            avg_err += j
        avg_err /= N

    print(f'Средняя ошибка в каждом N: {avg_err}')
    print(f'Средняя ошибка\t\t : {avg_err / K}')