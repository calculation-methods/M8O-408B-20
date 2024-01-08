from methods import *
from graphics import *


if __name__ == "__main__":

    print_parameters()

    lieb_u, lieb_iter = liebmann_method()
    lieb_err = solve_error(lieb_u)
    plot_liebmann(lieb_u, lieb_iter, lieb_err)

    seid_u, seid_iter = seidel_method()
    seid_err = solve_error(seid_u)
    plot_seidel(seid_u, seid_iter, seid_err)

    upper_relax_u, upper_relax_iter = upper_relaxation()
    upper_relax_err = solve_error(upper_relax_u)
    plot_upper_relaxation(upper_relax_u, upper_relax_iter, upper_relax_err)

    plot_analytical_solution(tt*hx)
    plt.show()
