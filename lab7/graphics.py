import matplotlib.pyplot as plt
from parameters import *

fig, ax = plt.subplots(1, 2, figsize=[12, 6])
fig.suptitle('lab 7 (var #10)', fontsize=15)

ax[0].set_title("methods")
ax[1].set_title("error")
ax[0].grid()
ax[1].grid()


def plot_analytical_solution(t):
    ax[0].plot(X, analytical_solution(X, t), label="analytical solution", color="blue")
    ax[0].legend()
    ax[1].legend()


def plot_liebmann(liebmann, liebmann_iters, liebmann_error):
    ax[0].plot(X, liebmann[:, tt], label=f"liebmann method (iters = {liebmann_iters})", color="red", linewidth=2)
    ax[1].plot(liebmann_error, label="liebmann error", color="red", linewidth=2)


def plot_seidel(seidel, seidel_iters, seidel_error):
    ax[0].plot(X, seidel[:, tt], label=f"seidel method (iters = {seidel_iters})", color="green", linewidth=4)
    ax[1].plot(seidel_error, label="seidel error", color="green", linewidth=4)


def plot_upper_relaxation(upper_relaxation, upper_relaxation_iters, upper_relaxation_error):
    ax[0].plot(X, upper_relaxation[:, tt], label=f"upper relaxation method (iters = {upper_relaxation_iters})", color="orange", linewidth=2)
    ax[1].plot(upper_relaxation_error, label="upper relaxation error", color="orange", linewidth=2)
