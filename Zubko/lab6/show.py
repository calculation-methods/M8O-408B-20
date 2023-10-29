import numpy as np
from matplotlib import pyplot as plt

from functions import analytic_solution
from task import count_t


def show_result(t_axis, x_axis, u1, u2):
    fig, ax = plt.subplots(2)
    fig.suptitle("Сравнение численных решений ДУ с аналитическим")
    fig.set_figheight(15)
    fig.set_figwidth(16)
    time = 0
    for i in range(2):
        ax[i].plot(x_axis, u1[time, :], label="Explicit scheme")
        ax[i].plot(x_axis, u2[time, :], label="Implicit scheme")
        ax[i].plot(
            x_axis,
            [analytic_solution(x, t_axis[time]) for x in x_axis],
            label="Analytic",
        )
        ax[i].grid(True)
        ax[i].set_xlabel("x")
        ax[i].set_ylabel("u")
        ax[i].set_title(f"Решения при t = {time / count_t}")
        time += count_t - 1

    plt.legend(bbox_to_anchor=(1.05, 2), loc="upper right", borderaxespad=0)
    plt.show()

    fig = plt.figure(num=1, figsize=(19, 12), clear=True)
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    fig.suptitle("Аналитическое решение")
    xgrid, tgrid = np.meshgrid(x_axis, t_axis)
    ax.plot_surface(xgrid, tgrid, analytic_solution(xgrid, tgrid))
    ax.set(xlabel="x", ylabel="t", zlabel="u")
    fig.tight_layout()
    plt.show()


def show_inaccuracy(t_axis, x_axis, u):
    inaccuracy = np.zeros(count_t)
    for i in range(count_t):
        inaccuracy[i] = np.max(
            np.abs(u[i] - np.array([analytic_solution(x, t_axis[i]) for x in x_axis]))
        )

    plt.figure(figsize=(14, 8))
    plt.plot(t_axis[1:], inaccuracy[1:], "violet", label="Ошибка")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper right", borderaxespad=0.0)
    plt.title("График изменения ошибки во времени")
    plt.xlabel("t")
    plt.ylabel("error")
    plt.grid(True)
    plt.show()
