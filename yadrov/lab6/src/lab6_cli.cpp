#include <iostream>
#include <cmath>
#include <iomanip>

#include ".././include/partial_differential/hyperbolic_pde.hpp"

using namespace HyperbolicPDE;

const double PI = 2 * std::acos(0.0);

int main() {
	PDE<double> pde(1, 0, 0, 0, [](double, double){return 0.;});
	pde.SetBoundaries(
		[](double x){return sin(PI * x);}, 
		[](double x){return PI * cos(PI * x);}, 
		[](double x){return -PI * PI * sin(PI * x);}, 
		[](double){return 0.;}, 
		0, 1,
		0, 1, [](double){return 0.;}, 
		0, 1, [](double){return 0.;});
	pde.SetSolution([&pde](double x, double t){return cos(PI * pde.a * t) * sin(PI * x);});

	int n = 10;
	double t_end = 1, sigma = 1;
	// std::cin >> t_end >> n >> sigma;

	std::vector<double> t, x;
	grid_t<double> u_expl, u_impl;
	std::cout << std::setprecision(12) << std::fixed;
	std::tie(x, t, u_expl) = ExplicitSolver<double>(pde, t_end, n, sigma, ApproxType::Taylor, ApproxType::Linear);
	u_impl = ImplicitSolver<double>(pde, x, t, t_end, ApproxType::Taylor, ApproxType::Linear);

	std::vector<double> error_impl(t.size()), error_expl(t.size());
    for (size_t k = 0; k < t.size(); ++k) {
        double max_error_impl = 0, max_error_expl = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            max_error_impl = std::max(max_error_impl, std::abs(pde.solution(x[i], t[k]) - u_impl[k][i]));
            max_error_expl = std::max(max_error_expl, std::abs(pde.solution(x[i], t[k]) - u_expl[k][i]));
        }
        error_expl[k] = max_error_expl;
        error_impl[k] = max_error_impl;
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
	std::cout << "EOD\n\n";

	std::cout << "set term qt size 1000, 800\n"
	             "set multiplot layout 2,2 rowsfirst\n\n"
	             "set size ratio -1\nset key off\n\n"
	             "set pm3d depthorder border lc \"white\" lt 1 lw 0.5\nset palette rgbformulae 33,13,10\n";

	std::cout << "set xlabel \"t\"\nset ylabel \"x\"\nset zlabel \"u(x, t)\"\n"
	             "set title \"Explicit scheme\"\n"
	             "splot $explicit with pm3d\n\n";

	std::cout << "set xlabel \"t\"\nset ylabel \"x\"\nset zlabel \"u(x, t)\"\n"
	             "set title \"Implicit scheme\"\n"
	             "splot $implicit with pm3d\n\n";

	std::cout << "set xlabel \"t\"\nset ylabel \"x\"\nset zlabel \"u(x, t)\"\n"
	             "set title \"Analitical solution\"\n"
	             "splot $solution with pm3d\n\n";

	std::cout << "set size square\nset key on\nset title \"Error\"\nset xlabel \"t\"\nset ylabel \"error\"\n"
	             "plot $errors every :::0::0 with lines title \"explicit\", \\\n"
	             "$errors every :::1::1 with lines title \"implicit\"\n\n";

	std::cout << "unset multiplot\n";
}
