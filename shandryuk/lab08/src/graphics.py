import numpy as np
import matplotlib.pyplot as plt

# Загрузка данных из файла
data = np.loadtxt('results.txt')

# Размеры сетки
Nx = 50
Ny = 50
Nt = 100

# Создание массивов для данных
x_values = data[:, 0]
y_values = data[:, 1]
t_values = data[:, 2]
numerical_solution = data[:, 3]
analytical_solution = data[:, 4]
error = data[:, 5]

# Построение графиков зависимости точности в разные моменты времени
unique_t_values = np.unique(t_values)
plt.figure(figsize=(15, 5))

for t in unique_t_values:
    indices = np.where(t_values == t)
    plt.plot(x_values[indices], error[indices], label=f'Time = {t:.2f}')

plt.title('Dependence of Accuracy on Time')
plt.xlabel('x')
plt.ylabel('Error')
plt.legend()
plt.grid(True)
plt.show()

# Построение графика зависимости погрешности от сеточных параметров
hx_values = np.unique(x_values)
hy_values = np.unique(y_values)

error_matrix = error.reshape((len(hx_values), len(hy_values)))

hx, hy = np.meshgrid(hx_values, hy_values)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(hx, hy, error_matrix.T, cmap='viridis')

ax.set_title('Dependence of Error on Grid Parameters')
ax.set_xlabel('hx')
ax.set_ylabel('hy')
ax.set_zlabel('Error')

plt.show()
