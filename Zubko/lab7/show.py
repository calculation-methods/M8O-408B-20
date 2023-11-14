import numpy as np
from matplotlib import pyplot as plt

from functions import analytic_solution
from task import count_x, count_y


def show_result(y_axis, x_axis, u1, u2, u3):
    fig, ax = plt.subplots(2)
    fig.suptitle("Сравнение численных решений ДУ с аналитическим")
    fig.set_figheight(15)
    fig.set_figwidth(16)
    y = 0
    for i in range(2):
        ax[i].plot(x_axis, u1[:, y], label="Liebmann method")
        ax[i].plot(x_axis, u2[:, y], label="Seidel method")
        ax[i].plot(x_axis, u3[:, y], label="Successive over-relaxation")
        ax[i].plot(
            x_axis, [analytic_solution(x, y_axis[y]) for x in x_axis], label="Analytic"
        )
        ax[i].grid(True)
        ax[i].set_xlabel("x")
        ax[i].set_ylabel("u")
        ax[i].set_title(f"Решения при y = {y / count_y}")
        y += count_y - 1

    plt.legend(bbox_to_anchor=(1.05, 2), loc="upper right", borderaxespad=0)
    plt.show()

    fig = plt.figure(num=1, figsize=(19, 12), clear=True)
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    fig.suptitle("Аналитическое решение")
    xgrid, ygrid = np.meshgrid(x_axis, y_axis)
    ax.plot_surface(xgrid, ygrid, analytic_solution(xgrid, ygrid))
    ax.set(xlabel="x", ylabel="y", zlabel="u")
    fig.tight_layout()
    plt.show()


def show_inaccuracy(y_axis, x_axis, u):
    inaccuracy = np.zeros(count_x + 1)
    for i in range(count_x + 1):
        inaccuracy[i] = np.max(
            np.abs(u[i] - np.array([analytic_solution(x_axis[i], y) for y in y_axis]))
        )

    plt.figure(figsize=(14, 8))
    plt.plot(x_axis[1:], inaccuracy[1:], "violet", label="Ошибка")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper right", borderaxespad=0.0)
    plt.title("График изменения ошибки")
    plt.xlabel("y")
    plt.ylabel("error")
    plt.grid(True)
    plt.show()
