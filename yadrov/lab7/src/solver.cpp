#include <iostream>
#include <vector>
#include <functional>
#include <cstring>
#include <cstdlib>
#include <tuple>
#include <iomanip>
#include <limits>

#include "../../include/partial_differential/elliptic_pde.hpp"
#include "presets.hpp"

using namespace EllipticPDE;

/*
Usage:
solver PRESET TYPE NX NY EPS RELAX
*/

int main(int argc, char const *argv[]) {
	if (argc < 7) {
		std::cerr << "Expected 6 arguments\n";
		exit(1);
	}

	std::map<int, PDE<double>> presets;
	GeneratePresets(presets);

	int preset, type, nx, ny;
	double eps, relax;

	char *end;
	preset = std::strtol(argv[1], &end, 10);
	if (end == argv[1] || presets.find(preset) == presets.end()) {
		std::cerr << "Invalid value for preset number " << preset << '\n';
		exit(1);
	}
	if (std::strcmp(argv[2], "iteration") == 0) {
		type = 0;
	} else if (std::strcmp(argv[2], "seidel") == 0) {
		type = 1;
	} else {
		std::cerr << "Invalid solver type\n";
		exit(1);
	}
	nx = std::strtol(argv[3], &end, 10);
	if (end == argv[3] || nx < 0) {
		std::cerr << "Invalid value for Nx\n";
		exit(1);
	}
	ny = std::strtol(argv[4], &end, 10);
	if (end == argv[4] || ny < 0) {
		std::cerr << "Invalid value for Ny\n";
		exit(1);
	}
	eps = std::strtod(argv[5], &end);
	if (end == argv[5] || eps < 0) {
		std::cerr << "Invalid value for eps\n";
		exit(1);
	}
	relax = std::strtod(argv[6], &end);
	if (end == argv[6] || relax < 0 || relax > 2) {
		std::cerr << "Invalid value for relax coefficient\n";
		exit(1);
	}

	std::vector<double> x, y;
	grid_t<double> u;
	int iter_count;
	PDE<double>& pde = presets[preset];

	if (type == 0) {
		std::tie(x, y, u, iter_count) = IterationSolver(pde, nx, ny, eps, relax);
	} else {
		std::tie(x, y, u, iter_count) = SeidelSolver(pde, nx, ny, eps, relax);
	}

	std::cout << iter_count << '\n';
	std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << std::fixed;
	for (int i = 0; i <= nx; ++i) {
		std::cout << x[i] << ' ';
	}
	std::cout << '\n';

	for (int i = 0; i <= ny; ++i) {
		std::cout << y[i] << ' ';
	}
	std::cout << '\n';

	for (int i = 0; i <= nx; ++i) {
		for (int j = 0; j <= ny; ++j) {
			std::cout << u[i][j] << ' ';
		}
	}
	std::cout << '\n';

	// for (int i = 0; i <= nx; ++i) {
	// 	for (int j = 0; j <= ny; ++j) {
	// 		std::cout << pde.solution(x[i], y[j]) << ' ';
	// 	}
	// }
	// std::cout << '\n';

	for (int i = 0; i <= nx; ++i) {
		for (int j = 0; j <= ny; ++j) {
			std::cout << pde.solution(x[i], y[j]) - u[i][j] << ' ';
		}
	}
}