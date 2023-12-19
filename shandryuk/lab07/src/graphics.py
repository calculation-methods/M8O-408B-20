import numpy as np
import matplotlib.pyplot as plt

def read_data(file_path):
    data = np.loadtxt(file_path, skiprows=1)
    x_values = data[:, 0]
    y_values = data[:, 1]
    numerical_solution = data[:, 2]
    analytical_solution = data[:, 3]
    error = data[:, 4]
    return x_values, y_values, numerical_solution, analytical_solution, error

def plot_error(error, title):
    plt.figure(figsize=(8, 6))
    plt.plot(error, label='Error')
    plt.title(title)
    plt.xlabel('Data Point')
    plt.ylabel('Error')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_error_dependency(h_values, error, title, xlabel):
    plt.figure(figsize=(8, 6))
    plt.plot(h_values, error, marker='o')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Max Error')
    plt.grid(True)
    plt.show()

# Чтение данных из файлов
x_values_libman, y_values_libman, _, _, error_libman = read_data('result_libman.txt')
x_values_zeidel, y_values_zeidel, _, _, error_zeidel = read_data('result_zeidel.txt')
x_values_prost_iter, y_values_prost_iter, _, _, error_prost_iter = read_data('result_prost_iter.txt')

# Построение графиков погрешности
plot_error(error_libman, 'Error Plot - Libman Method')
plot_error(error_zeidel, 'Error Plot - Zeidel Method')
plot_error(error_prost_iter, 'Error Plot - Prost Iter Method')

# Зависимость погрешности от hx и hy
hx_values_libman = np.unique(x_values_libman)
hy_values_libman = np.unique(y_values_libman)
hx_values_zeidel = np.unique(x_values_zeidel)
hy_values_zeidel = np.unique(y_values_zeidel)
hx_values_prost_iter = np.unique(x_values_prost_iter)
hy_values_prost_iter = np.unique(y_values_prost_iter)

plot_error_dependency(hx_values_libman, error_libman, 'Error vs hx - Libman Method', 'hx')
plot_error_dependency(hy_values_libman, error_libman, 'Error vs hy - Libman Method', 'hy')

plot_error_dependency(hx_values_zeidel, error_zeidel, 'Error vs hx - Zeidel Method', 'hx')
plot_error_dependency(hy_values_zeidel, error_zeidel, 'Error vs hy - Zeidel Method', 'hy')

plot_error_dependency(hx_values_prost_iter, error_prost_iter, 'Error vs hx - Prost Iter Method', 'hx')
plot_error_dependency(hy_values_prost_iter, error_prost_iter, 'Error vs hy - Prost Iter Method', 'hy')
