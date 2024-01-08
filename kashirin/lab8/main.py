from methods import *
from graphics import *

if __name__ == "__main__":

    print_parameters()

    u_mdn = mdn()
    err_mdn = solve_error(u_mdn)
    plot_mdn(u_mdn, err_mdn)

    u_mds = mds()
    err_mds = solve_error(u_mds)
    plot_mds(u_mds, err_mds)

    plot_analytical_solution()

    plt.show()

