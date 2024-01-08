import matplotlib.pyplot as plt
import numpy as np

# Пример №1
kernel = lambda x, s: np.exp(s - x)
right_side = lambda x: np.exp(-x)
exact_solution = lambda x: 1 + x * 0 
start = 0
end = 1
step_size = 0.05

def quad_trapezoid_method(start, end, step_size, kernel, right_side):
    num_intervals = int((end - start) / step_size)
    x_values = np.linspace(start, end, num_intervals + 1)
    y_values = np.empty_like(x_values)
    y_values[0] = right_side(x_values[0])
    kernel_ii = kernel(x_values[0], x_values[0])
    
    for i in range(1, num_intervals + 1):
        c_i = right_side(x_values[i]) + sum(step_size * kernel(x_values[i], x_values[j]) * y_values[j] * y_values[j] for j in range(0, i))
        y_values[i] = (1 - np.sqrt(1 - 2 * step_size * kernel_ii * c_i)) / (step_size * kernel_ii)
    
    return x_values, y_values

def runge_romberg_error_estimation(method, order, step_size):
    x_h, y_h = method(start, end, step_size, kernel, right_side)
    x_2h, y_2h = method(start, end, 2 * step_size, kernel, right_side)
    num_intervals = len(x_2h)
    error = max(abs((y_h[2 * i] - y_2h[i])) / (2 ** order - 1) for i in range(num_intervals))
    print("Рунге-Ромберг =", error)

figsize_result = (8, 6)
figsize_error = (8, 4)

# Результат метода
x_result, y_result = quad_trapezoid_method(start, end, step_size, kernel, right_side)
plt.figure(figsize=figsize_result)
plt.plot(x_result, y_result, 'o', label='Результат', color='blue')
plt.plot(x_result, exact_solution(x_result), '-', lw=2, label='Решение', color='green')
plt.grid(True)
plt.title('Метод квадратур')
plt.show()

# График ошибки
plt.figure(figsize=figsize_error)
plt.grid(True)
plt.title('График ошибки')
plt.plot(x_result, abs(y_result - exact_solution(x_result)), '-', color='red')
plt.show()

max_error = max(abs(y_result - exact_solution(x_result)))
print("Ошибка = ", max_error)

runge_romberg_error_estimation(quad_trapezoid_method, 2, step_size)






# Пример №2
kernel = lambda x, s: 1 + x*s*0
right_side = lambda x: np.sin(x) - x/2 + np.sin(2*x)/4
exact_solution = lambda x: np.sin(x)
start = 0
end = np.pi/2
step_size = np.pi/2 / 14

# Результат метода
x_result, y_result = quad_trapezoid_method(start, end, step_size, kernel, right_side)
plt.figure(figsize=figsize_result)
plt.plot(x_result, y_result, 'o', label='Результат', color='blue')
plt.plot(x_result, exact_solution(x_result), '-', lw=2, label='Решение', color='green')
plt.grid(True)
plt.title('Метод квадратур')
plt.show()

# График ошибки
plt.figure(figsize=figsize_error)
plt.grid(True)
plt.title('График ошибки')
plt.plot(x_result, abs(y_result - exact_solution(x_result)), '-', color='red')
plt.show()

max_error = max(abs(y_result - exact_solution(x_result)))
print("Ошибка = ", max_error)

runge_romberg_error_estimation(quad_trapezoid_method, 2, step_size)



# Пример №3
kernel = lambda x, s, y: np.cos(np.pi*x) * np.sin(np.pi*s) * y**3 / 5
kernel_derivative = lambda x, s, y: np.cos(np.pi*x) * np.sin(np.pi*s) * y*y * 3 / 5
right_hand_side = lambda x: np.sin(np.pi*x)
true_solution = lambda x: np.sin(np.pi*x) + np.cos(np.pi*x)*(20 - np.sqrt(391)) / 3

# Параметры задачи
a = 0
b = 1
h = 0.05
m = 2

def runge_romberg_2(method, p, h):
    x_h, y_h = method(kernel, kernel_derivative, right_hand_side, a, b, h, m)
    x_2h, y_2h = method(kernel, kernel_derivative, right_hand_side, a, b, 2*h, m)
    N = len(x_2h)
    R = max(abs((y_h[2*i] - y_2h[i])) / (2**p - 1) for i in range(N))
    print("Рунге-Ромберг = ", R)

def newton_kantorovich(K, K_y, f, a, b, h, m):
    k = 0
    N = int((b - a) / h)
    x = np.linspace(a, b, N+1)
    I = np.eye(N+1)
    Y = np.ones_like(x)
    F = np.empty_like(Y)
    A = np.empty((N+1, N+1))
    
    while k < m: 
        for i in range(N+1):
            F[i] = f(x[i]) + sum(h * (K(x[i], x[j], Y[j]) - K_y(x[i], x[j], Y[j]) * Y[j]) for j in range(1, N)) + h/2 * (
                        K(x[i], x[0], Y[0]) + K(x[i], x[N], Y[N]) - K_y(x[i], x[0], Y[0]) * Y[0] - K_y(x[i], x[N], Y[N]) * Y[N])
        for i in range(N+1):
            A[i, 0] = K_y(x[i], x[0], Y[0]) * h/2
            A[i, N] = K_y(x[i], x[N], Y[N]) * h/2
            for j in range(1, N):
                A[i, j] = K_y(x[i], x[j], Y[j]) * h
        Y = np.linalg.solve(I - A, F)
        k += 1
    
    return x, Y

x_result, y_result = newton_kantorovich(kernel, kernel_derivative, right_hand_side, a, b, h, m)
plt.plot(x_result, y_result, 'o', label='Результат', color='blue')
plt.plot(x_result, true_solution(x_result), '-', lw=2, label='Решение', color='green')
plt.grid(True)
plt.title('Ньютон-Канторович')
plt.show()

plt.grid(True)
plt.title('График ошибки')
plt.plot(x_result, abs(y_result - true_solution(x_result)), '-', color='red')
plt.show()
print("Ошибка = ", max(abs(y_result - true_solution(x_result))))

runge_romberg_2(newton_kantorovich, 2, h)