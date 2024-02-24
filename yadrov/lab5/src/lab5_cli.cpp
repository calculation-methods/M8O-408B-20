#include <iostream>
#include <cmath>
#include <iomanip>

#include "../../include/partial_differential/parabolic_pde.hpp"

using namespace ParabolicPDE;

const double PI = 2 * std::acos(0.0);

int main() {
	PDE<double> pde{
		1., 1., -1., [](double, double){ return 0; }
	};
	pde.SetBoundaries(
		[](double x){ return sin(x); },
		0., PI,
		1., 1., [&pde](double t){ return exp((pde.c - pde.a) * t) * (cos(pde.b * t) + sin(pde.b * t)); },
		1., 1., [&pde](double t){ return -exp((pde.c - pde.a) * t) * (cos(pde.b * t) + sin(pde.b * t)); });
	pde.SetSolution(
		[&pde](double x, double t){ return exp((pde.c - pde.a) * t) * sin(x + pde.b * t); });

	int n;
	double t_end, sigma, theta = 0.5;
	std::cin >> t_end >> n >> sigma;

	std::vector<double> t, x;
	grid_t<double> u_expl, u_impl, u_comb;
	std::cout << std::setprecision(12) << std::fixed;
	std::tie(x, t, u_expl) = ExplicitSolver<double>(pde, t_end, n, sigma, ApproxType::Taylor);
	u_impl = ImplicitSolver<double>(pde, x, t, t_end, ApproxType::Taylor);
	u_comb = CombinedSolver<double>(pde, x, t, t_end, ApproxType::Taylor, theta);

	std::vector<double> error_impl(t.size()), error_expl(t.size()), error_comb(t.size());
    for (size_t k = 0; k < t.size(); ++k) {
        double max_error_impl = 0, max_error_expl = 0, max_error_comb = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            max_error_impl = std::max(max_error_impl, std::abs(pde.solution(x[i], t[k]) - u_impl[k][i]));
            max_error_expl = std::max(max_error_expl, std::abs(pde.solution(x[i], t[k]) - u_expl[k][i]));
            max_error_comb = std::max(max_error_comb, std::abs(pde.solution(x[i], t[k]) - u_comb[k][i]));
        }
        error_expl[k] = max_error_expl;
        error_impl[k] = max_error_impl;
        error_comb[k] = max_error_comb;
    }

	std::cout << "$explicit << EOD\n";
	for (size_t k = 0; k < t.size(); ++k) {
		for (size_t i = 0; i < x.size(); ++i) {
			std::cout << t[k] << ' ' << x[i] << ' ' << u_expl[k][i] << '\n';
		}
		std::cout << '\n';
	}
	std::cout << "EOD\n";

	std::cout << "$implicit << EOD\n";
	for (size_t k = 0; k < t.size(); ++k) {
		for (size_t i = 0; i < x.size(); ++i) {
			std::cout << t[k] << ' ' << x[i] << ' ' << u_impl[k][i] << '\n';
		}
		std::cout << '\n';
	}
	std::cout << "EOD\n";

	std::cout << "$combined << EOD\n";
	for (size_t k = 0; k < t.size(); ++k) {
		for (size_t i = 0; i < x.size(); ++i) {
			std::cout << t[k] << ' ' << x[i] << ' ' << u_comb[k][i] << '\n';
		}
		std::cout << '\n';
	}
	std::cout << "EOD\n";

	std::cout << "$solution << EOD\n";
	for (size_t k = 0; k < t.size(); ++k) {
		for (size_t i = 0; i < x.size(); ++i) {
			std::cout << t[k] << ' ' << x[i] << ' ' << pde.solution(x[i], t[k]) << '\n';
		}
		std::cout << '\n';
	}
	std::cout << "EOD\n";

	std::cout << "$errors << EOD\n";
	for (size_t k = 0; k < t.size(); ++k) {
		std::cout << t[k] << ' ' << error_expl[k] << '\n';
	}
	std::cout << '\n';
	for (size_t k = 0; k < t.size(); ++k) {
		std::cout << t[k] << ' ' << error_impl[k] << '\n';
	}
	std::cout << '\n';
	for (size_t k = 0; k < t.size(); ++k) {
		std::cout << t[k] << ' ' << error_comb[k] << '\n';
	}
	std::cout << "EOD\n\n";

	std::cout << "set term qt size 1500, 1000\n"
	             "set multiplot layout 2,3 rowsfirst\n\n"
	             "set size ratio -1\nset key off\n\n"
	             "set pm3d depthorder border lc \"white\" lt 1 lw 0.5\nset palette rgbformulae 33,13,10\n";

	std::cout << "set xlabel \"t\"\nset ylabel \"x\"\nset zlabel \"u(x, t)\"\n"
	             "set title \"Explicit scheme\"\n"
	             "splot $explicit with pm3d\n\n";

	std::cout << "set xlabel \"t\"\nset ylabel \"x\"\nset zlabel \"u(x, t)\"\n"
	             "set title \"Implicit scheme\"\n"
	             "splot $implicit with pm3d\n\n";

	std::cout << "set xlabel \"t\"\nset ylabel \"x\"\nset zlabel \"u(x, t)\"\n"
	             "set title \"Combined scheme\"\n"
	             "splot $combined with pm3d\n\n";

	std::cout << "set xlabel \"t\"\nset ylabel \"x\"\nset zlabel \"u(x, t)\"\n"
	             "set title \"Analitical solution\"\n"
	             "splot $solution with pm3d\n\n";

	std::cout << "set size square\nset key on\nset title \"Error\"\nset xlabel \"t\"\nset ylabel \"error\"\n"
	             "plot $errors every :::0::0 with lines title \"explicit\", \\\n"
	             "$errors every :::1::1 with lines title \"implicit\", \\\n"
	             "$errors every :::2::2 with lines title \"combined\"\n\n";

	std::cout << "unset multiplot\n";
}
