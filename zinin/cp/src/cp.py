import numpy as np
import matplotlib.pyplot as plt

def f(x, y, z):
    return -0.04 * x + 1.4 * y * z

def g(x, y, z):
    return 0.04 * x - 1.4 * y * z - 3 * y**2

def h(x, y, z):
    return -0.1 * z + 1 + x

def newton_method(F, J, y0, tol=1e-10, max_iter=1000):
    y = np.array(y0)
    for _ in range(max_iter):
        F_val = F(y)
        J_val = J(y)
        y_new = y - np.linalg.solve(J_val, F_val)
        if np.linalg.norm(y_new - y) < tol:
            return y_new
        y = y_new
    raise RuntimeError("Метод Ньютона не сошелся")

def implicit_runge_kutta(F, y0, t0, tf, h):
    def G(y, y_prev, t, h):
        return y - y_prev - h/2 * (F(t, y) + F(t - h, y_prev))

    def J_G(y, t, h):
        return np.eye(len(y)) - h/2 * J_F(t, y)

    t_values = [t0]
    y_values = [y0]

    y = y0
    t = t0
    while t < tf:
        t += h
        y = newton_method(lambda y_new: G(y_new, y, t, h), lambda y_new: J_G(y_new, t, h), y)
        t_values.append(t)
        y_values.append(y)

    return t_values, y_values

def F(t, y):
    x, y, z = y
    return np.array([f(x, y, z), g(x, y, z), h(x, y, z)])

def J_F(t, y):
    x, y, z = y
    # Якобиан для трехмерного случая
    return np.array([[-0.04, 1.4*z, 1.4*y],
                     [0.04, -1.4*z - 6*y, -1.4*y],
                     [1, 0, -0.1]])

# Начальные условия для x, y, z
y0 = np.array([1, 0.1, 0])  # Трехмерный вектор

t0, tf = 0, 10
hh = 0.01

# Вычисление численного решения
t_values, y_values = implicit_runge_kutta(F, y0, t0, tf, hh)

# Преобразование списка значений y в массив NumPy
y_values = np.array(y_values)

# Визуализация для всех трех компонент
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.plot(t_values, y_values[:, 0], label='x(t)')
plt.xlabel('Time')
plt.ylabel('x(t)')
plt.legend()

plt.subplot(1, 3, 2)
plt.plot(t_values, y_values[:, 1], label='y(t)')
plt.xlabel('Time')
plt.ylabel('y(t)')
plt.legend()

plt.subplot(1, 3, 3)
plt.plot(t_values, y_values[:, 2], label='z(t)')
plt.xlabel('Time')
plt.ylabel('z(t)')
plt.legend()

plt.tight_layout()
plt.show()
