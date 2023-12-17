import numpy as np
import matplotlib.pyplot as plt

from solver import ParabolicSolver
from utils import *
from config import *

def PlotUX(dict_, nx, ny, T, K, time, method_name="Метод дробных шагов"):
    hx = 1 / nx
    hy = 1 / ny
    tau = T / K

    x = np.arange(0, 1 + hx, hx)
    y = np.arange(0, 1 + hy, hy)
    t = np.arange(0, T + tau, tau)

    z1 = []
    for i in range(len(x)):
        z1.append(approximation(x, 10, time, dict_))

    z2 = solving(x, y[10], t[time])

    # Plot for analytic solution
    fig2 = plt.figure()
    plt.title(f'U(x) {method_name}')
    plt.plot(x, z2, color='b', label='аналитический')
    plt.plot(x, z1[time], color='r', label='численный')
    plt.legend(loc='best')
    plt.ylabel('U')
    plt.xlabel('x')
    plt.savefig('u(x).png')
    plt.show()
    

def plotDependenceError(Nx=4, Ny=4, K=40, Method=0, save_file="dependence_error_plot.png"):
    NxFix = Nx
    NyFix = Ny
    if Method != 0:
        method = True
        stri = 'Зависимость погрешности от шага (метода переменных направлений)'
    else:
        method = False
        stri = 'Зависимость погрешности от шага (метода дробных шагов)'
    plt.figure(figsize=(8, 8))
    Time = TimeList[2-1]  # random.randint(0, K - 1)
    schema = Schema(T=Time, order2nd=method)
    h, eps = [], []
    for j in range(10):
        h.append([])
        eps.append([])
        X, Y, T, Z = schema(Nx, Ny, K)
        Nx += 1
        Ny += 1
        h[-1].append(schema.hx)
        eps[-1].append(Error(X, Y, T[Time], Z[Time]))
    Nx = NxFix
    Ny = NyFix
    plt.plot(np.array(h), np.array(eps))
    plt.xlabel('$h$')
    plt.ylabel('погрешность')
    plt.title(stri)
    plt.grid()

    # Save the figure
    plt.savefig(save_file)

    
def PlotError(save_file="error_surface_plot.png"):
    schema = Schema(T=5, order2nd=True)
    h = []
    tau = []
    eps = []

    for i in range(20):
        h.append([])
        tau.append([])
        eps.append([])
        for j in range(40):
            N = i + 5
            K = j + 40
            X, Y, T, Z = schema(N, N, K)
            h[-1].append(schema.hx)
            tau[-1].append(schema.tau)
            eps[-1].append(FullError(X, Y, T, Z))

    fig = plt.figure(num=1, figsize=(19, 12), clear=True)
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.plot_surface(np.array(h), np.array(tau), np.array(eps))
    ax.set(xlabel='$h_x$', ylabel='$t$', zlabel='$\epsilon$', title='Погрешность метода переменных направлений')
    fig.tight_layout()

    plt.savefig(save_file)
    
def PlotTimeError(Method=0, save_file="error_time_plot.png"):
	


	if Method != 0:
		method = True
		stri = 'Зависимость погрешности от времени (метода переменных направлений)'
	else:
		method = False
		stri = 'Зависимость погрешности от времени (метода дробных шагов)'
        
	schema = Schema(T=5, order2nd=method)
	h = []
	tau = []
	eps = []

	for i in range(20):
		h.append([])
		tau.append([])
		eps.append([])
		for j in range(40):
			N = i + 5
			K = j + 40
			X, Y, T, Z = schema(N, N, K)
			h[-1].append(schema.hx)
			tau[-1].append(schema.tau)
			eps[-1].append(FullError(X, Y, T, Z))

	t = range(0, 100, 5)
	err = [max(m) for m in eps]
	fig = plt.figure(figsize=(9, 9))
	plt.title(stri)
	plt.plot(t, err)
	plt.legend(loc='best')
	plt.ylabel('погрешность')
	plt.xlabel('t')
	plt.show()
	plt.savefig(save_file)
    

if __name__ == '__main__':

	solver = ParabolicSolver(params, equation_type)
	solved = solver.solve(nx, ny, T, K)
	ans = {
		'numerical': solved,
		'analytic': solver.analyticSolve(nx, ny, T, K)
	}
	plotDependenceError(Method=1)
	plotDependenceError(Method=0)
	PlotUX(ans, nx, ny, T, K, 20)
	PlotTimeError(Method=0)
	PlotTimeError(Method=1)

	
