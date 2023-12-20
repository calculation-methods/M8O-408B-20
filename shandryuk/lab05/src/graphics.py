import numpy as np
import matplotlib.pyplot as plt

# Считывание данных из файлов
def read_data(file_path):
    data = np.loadtxt(file_path)
    x = data[:, 0]
    t = data[:, 1]
    u_numerical = data[:, 2]
    u_analytical = data[:, 3]
    return x, t, u_numerical, u_analytical

# Построение графика
def plot_graph(x, t, u_numerical, u_analytical, title):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, t, u_numerical, label='Numerical Solution', color='r', marker='o')
    ax.scatter(x, t, u_analytical, label='Analytical Solution', color='b', marker='^')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u')
    plt.title(title)
    plt.legend()
    plt.show()

# Чтение и построение графика для каждого файла
files = ["results_explicit.txt", "results_implicit.txt", "results_cn.txt"]

for file_path in files:
    x, t, u_numerical, u_analytical = read_data(file_path)
    title = f'Numerical vs Analytical Solution\n({file_path})'
    plot_graph(x, t, u_numerical, u_analytical, title)
