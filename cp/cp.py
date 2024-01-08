import numpy as np
import math
import plotly.offline as offline
from plotly.graph_objs import *

def solution(x):
  return math.exp(x**2 + x)

def f_volt(x):
  return math.exp(x**2)

def K_volt(x, s):
  return math.exp(x**2 - s**2)

def linear_volterra_2(a, b, n, h, x, f, K):
    phi = [f(a)]
    
    for i in range(1, n):
        s = 0
        j = 0
        while (j < i):
            s += K(x[i], x[j]) * phi[j]
            j += 1
        phi.append((1 / (1 - h/2 * K(x[i], x[i]))) * (f(x[i]) + h/2 * K(x[i], x[0]) * phi[0] + h * s))
    return x, phi

a = 0.0
b = 1.0
n = 52

h = (b - a) / (n - 1)
    
x = list(np.arange(a, b + h, h))
print("Пример 1: ")
sol_exact = [solution(i) for i in x]
print("Точное решение: ")
print(sol_exact)

x, sol_sq1 = linear_volterra_2(a, b, n, h, x, f_volt, K_volt)
print("Метод квадратур: ")
print(sol_sq1)

def obs_error(sol1, sol2, x_list):
  res = math.sqrt(sum((sol1[i] - sol2[i])**2 for i in range(len(x_list))))
  return res

print("Погрешность: ")
print(obs_error(sol_sq1, sol_exact, x))

scatter_exact = Scatter(x = x, y = sol_exact, name = 'Точное решение', mode = 'markers + lines', showlegend = True)
scatter_sq = Scatter(x = x, y = sol_sq1, name = 'Метод квадратур', mode = 'markers + lines', showlegend = True)
data = [scatter_exact, scatter_sq]
layout = Layout(xaxis = dict(title = 'x'), yaxis = dict(title = 'y'))
fig = Figure(data = data, layout = layout)
offline.plot(fig)

def solution2(x):
  return 0.5 * x * (math.exp(-x) + math.exp(x))

def f_volt2(x):
  return x

def K_volt2(x, s):
  return (4 * math.sin(x - s) - x + s)

a = 0.0
b = 1.0
n = 32

h = (b - a) / (n - 1)
    
x = list(np.arange(a, b + h, h))
print("Пример 2: ")
sol_exact2 = [solution2(i) for i in x]
print("Точное решение: ")
print(sol_exact2)

x, sol_sq2 = linear_volterra_2(a, b, n, h, x, f_volt2, K_volt2)
print("Метод квадратур: ")
print(sol_sq2)

print("Погрешность: ")
print(obs_error(sol_sq2, sol_exact2, x))

scatter_exact = Scatter(x = x, y = sol_exact2, name = 'Точное решение', mode = 'markers + lines', showlegend = True)
scatter_sq = Scatter(x = x, y = sol_sq2, name = 'Метод квадратур', mode = 'markers + lines', showlegend = True)
data = [scatter_exact, scatter_sq]
layout = Layout(xaxis = dict(title = 'x'), yaxis = dict(title = 'y'))
fig = Figure(data = data, layout = layout)
offline.plot(fig)