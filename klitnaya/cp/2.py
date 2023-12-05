from math import sin, cos, exp, pi, log, tan, asin, acos, atan
import math
import numpy as np

# Функция для считывания пользовательской функции
def read_function(prompt):
    while True:
        try:
            # Ввод пользовательской функции
            function_str = input(prompt)
            # Вычисление функции с использованием eval()
            function = eval('lambda x, y: ' + function_str)
            break
        except:
            print("Ошибка! Пожалуйста, введите корректную функцию.")
    return function

def read_functionx(prompt):
    while True:
        try:
            # Ввод пользовательской функции
            function_str = input(prompt)
            # Вычисление функции с использованием eval()
            function = eval('lambda x: ' + function_str)
            break
        except:
            print("Ошибка! Пожалуйста, введите корректную функцию.")
    return function
'''
# Задание функции ядра
def kernel(x, y):
    # TODO: Вставьте код для вычисления функции ядра
    return x*y

# Задание функции f
def f(x):
    # TODO: Вставьте код для вычисления функции f
    return 0
'''
while True:
    # Задание границы интегрирования
    print("введите границы интегрирования a и b: ")
    a = float(input())#0
    b = float(input())#1

    print("Введите количества узлов")

    # Задание количества узлов
    n = int(input())#100

    # Расчет шага интегрирования
    h = (b - a) / n

    # Считывание функции ядра
    kernel_prompt = "Введите функцию ядра (от переменных x и y): "
    kernel = read_function(kernel_prompt)

    # Считывание функции f
    f_prompt = "Введите функцию f (от переменной x): "
    f = read_functionx(f_prompt)

    # Вычисление матрицы A
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            x = a + i * h
            y = a + j * h
            A[i][j] = kernel(x, y)

    # Вычисление вектора B
    B = [0] * n
    for i in range(n):
        x = a + i * h
        B[i] = f(x)

    # Решение системы линейных уравнений
    x = [0] * n
    for i in range(n):
        for j in range(n):
            x[i] += A[i][j] * B[j]

    # Вычисление приближенного значения интеграла
    integral = 0
    for i in range(n):
        integral += x[i] * h

    # Вывод результата
    print("Приближенное значение интеграла:", integral)
    choice = input("Желаете продолжить? (y/n): ")
    if choice.lower() != "y":
        break

