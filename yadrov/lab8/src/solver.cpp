#include <iostream>
#include <vector>
#include <functional>
#include <tuple>
#include <iomanip>
#include <limits>
#include <cmath>

#include "../../include/partial_differential/parabolic2d_pde.hpp"

using namespace Parabolic2dPDE;

using std::sin, std::cos, std::exp, std::log, std::cosh, std::sinh;
const double pi = 2 * std::acos(0);

int main() {
	// вариант 10
	double mu = pi;
	PDE<double> pde{
		1, 0, 0, 0, [mu](double x, double y, double t){ return sin(x) * sin(y) * (mu * cos(mu * t) + 2 * sin(mu * t)); },
		0, 2 * pi, 0, 2 * pi,
		[](double, double){ return 0; },
		0, 1, [](double, double){ return 0; },
		1, 0, [mu](double y, double t){ return sin(y) * sin(mu * t); },
		0, 1, [](double, double){ return 0; },
		1, 0, [mu](double x, double t){ return sin(x) * sin(mu * t); },
		[mu](double x, double y, double t){ return sin(x) * sin(y) * sin(mu * t); }
	};

	int method, nx, ny, nt;
	double t_end;
	std::cin >> method >> nx >> ny >> nt >> t_end;

	std::vector<double> x, y, t;
	Tensor<3, double> u;

	if (method == 0) {
		std::tie(u, x, y, t) = AlternatingDirectionMethod(pde, nx, ny, nt, t_end);
	} else {
		std::tie(u, x, y, t) = FractionalStepMethod(pde, nx, ny, nt, t_end);
	}

	std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << std::fixed;
	for (int i = 0; i <= nx; ++i) {
		std::cout << x[i] << ' ';
	}
	std::cout << '\n';

	for (int i = 0; i <= ny; ++i) {
		std::cout << y[i] << ' ';
	}
	std::cout << '\n';

	for (int i = 0; i <= nt; ++i) {
		std::cout << t[i] << ' ';
	}
	std::cout << '\n';

	for (int k = 0; k <= nt; ++k) {
		for (int i = 0; i <= nx; ++i) {
			for (int j = 0; j <= ny; ++j) {
				std::cout << u[k][i][j] << ' ';
			}
		}
	}
	std::cout << '\n';

	for (int k = 0; k <= nt; ++k) {
		for (int i = 0; i <= nx; ++i) {
			for (int j = 0; j <= ny; ++j) {
				std::cout << pde.solution(x[i], y[j], t[k]) - u[k][i][j] << ' ';
			}
		}
	}
}