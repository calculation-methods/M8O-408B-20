import matplotlib.pyplot as plt
from parameters import *

fig, ax = plt.subplots(1, 2, figsize=[12, 9])
fig.suptitle('lab 8 (var #10)', fontsize=15)

ax[0].set_title("methods")
ax[1].set_title("error")
ax[0].grid()
ax[1].grid()


def plot_analytical_solution():
    ax[0].plot(y, analytical_solution(x[kx], y, t[kt]), label="analytical solution", color="blue")
    ax[0].legend()
    ax[1].legend()


def plot_mdn(u_mdn, u_mdn_error):
    ax[0].plot(y, u_mdn[kx, :, kt], label="mpn", color="red")
    ax[1].plot(x, u_mdn_error, label="mpn error", color="red")


def plot_mds(u_mds, u_mds_error):
    ax[0].plot(y, u_mds[kx, :, kt], label="mds", color="green")
    ax[1].plot(x, u_mds_error, label="mds error", color="green")
