import numpy as np
from matplotlib import pyplot as plt

from functions import analytic_solution
from task import count_x, count_y, hx, hy, lx, ly, t_count, t_max, tau


def show_result(u1, u2):
    x = np.arange(0, lx + hx, hx)
    y = np.arange(0, ly + hy, hy)
    t = np.arange(0, t_max + tau, tau)
    x_i = 1
    t_i = 2
    fig, ax = plt.subplots(2)
    fig.suptitle("Сравнение численных решений ДУ с аналитическим")
    fig.set_figheight(15)
    fig.set_figwidth(16)

    for i in range(2):
        ax[i].plot(
            y, analytic_solution(x[x_i], y, t[t_i]), color="red", label="Analytic"
        )
        ax[i].plot(y, u1[x_i, :, t_i], label="Схема переменных направлений")
        ax[i].plot(y, u2[x_i, :, t_i], label="Схема дробных шагов")
        ax[i].grid(True)
        ax[i].set_xlabel("y")
        ax[i].set_ylabel("u")
        ax[i].set_title(f"Решения при x= {x[x_i]}, t = {t[t_i]}")
        x_i = count_x
        t_i = t_count

    plt.legend(bbox_to_anchor=(1.05, 2), loc="upper right", borderaxespad=0)
    plt.show()


def show_inaccuracy(u1, u2):
    x = np.arange(0, lx + hx, hx)
    y = np.arange(0, ly + hy, hy)
    t = np.arange(0, t_max + tau, tau)
    x_i = 1
    t_i = 2
    plt.title(f"Погрешность по y при x= {x[x_i]}, t = {t[t_i]}")

    for i in range(2):
        inaccuracy = []
        if i == 0:
            u_tfix = u1[:, :, t_i]
        else:
            u_tfix = u2[:, :, t_i]

        for j in range(count_y + 1):
            a = analytic_solution(x, y[j], t[t_i]) - u_tfix[:, j]
            inaccuracy = np.append(inaccuracy, np.linalg.norm(a))
        plt.plot(y, inaccuracy)

    plt.xlabel("y")
    plt.ylabel("error")
    plt.xlim((0, ly))
    plt.grid(True)
    plt.show()
