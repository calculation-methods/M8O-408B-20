#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "../include/function/exponential_spline.hpp"

/*
-tense COUNT - number of iterations
-relax VALUE - relaxation parameter
-p VALUE - default tension parameter
-plot - gnuplot script
*/

int main(int argc, char* argv[]) {
	bool plot = false;
	bool default_p = false;
	double p_value = 0;
	double relax = 1;
	unsigned long long tense = 0;

	char *end;
	for (int i = 1; i < argc; ++i) {
		if (!plot && std::strcmp(argv[i], "-plot") == 0) {
			plot = true;

		} else if (!default_p && std::strcmp(argv[i], "-p") == 0) {
			default_p = true;
			if (++i == argc) {
				std::cerr << "Expected value after -p\n";
				return 1;
			}
			p_value = std::strtod(argv[i], &end);
			if (end == argv[i] || p_value <= 0 || p_value == HUGE_VAL) {
				std::cerr << "Invalid -p option value: " << argv[i] << '\n';
				return 1;
			}

		} else if (std::strcmp(argv[i], "-relax") == 0) {
			if (++i == argc) {
				std::cerr << "Expected value after -relax\n";
				return 1;
			}
			relax = std::strtod(argv[i], &end);
			if (end == argv[i] || relax < 0 || relax == HUGE_VAL) {
				std::cerr << "Invalid -relax option value: " << argv[i] << '\n';
				return 1;
			}

		} else if (std::strcmp(argv[i], "-tense") == 0) {
			if (++i == argc) {
				std::cerr << "Expected value after -tense\n";
				return 1;
			}
			tense = std::strtoull(argv[i], &end, 10);
			if (end == argv[i] || tense == ULLONG_MAX) {
				std::cerr << "Invalid -tense option value: " << argv[i] << '\n';
				return 1;
			}

		} else {
			std::cerr << "Unknown option: " << argv[i] << std::endl;
		}
	}

	std::vector<double> x, y, p;
	double a, b, c;

	while (std::cin >> a >> b) {
		x.push_back(a);
		y.push_back(b);
		if (!default_p && std::cin >> c) {
			p.push_back(c);
		}
	}

	if (default_p) {
		p.assign(x.size() - 1, p_value);
	}
	
	ExponentialSpline<double> spline(x, y, p);

	for (unsigned long long i = 0; i < tense; ++i) {
		spline.Tense(relax);
	}

	if (!plot) {
		std::cout << spline;

	} else {
		std::cout << "set term pdfcairo\nset size ratio -1\nset key off\n";
		const auto [xmin, xmax] = std::minmax_element(spline.x.begin(), spline.x.end());
		const auto [ymin, ymax] = std::minmax_element(spline.y.begin(), spline.y.end());
		double dx = (*xmax - *xmin) * 0.1, dy = (*ymax - *ymin) * 0.5;
		std::cout << "set xrange [" << (*xmin-dx) << ':' << (*xmax+dx) << "]\nset yrange [" << (*ymin-dy) << ':' << (*ymax+dy) << "]\n";
		std::cout << "set style line 1 lw 1.5 lt 1 pt 7 ps 0.5\n";
		std::cout << "$points << EOD\n";
		for (size_t i = 0; i < x.size(); ++i) {
			std::cout << x[i] << ' ' << y[i] << '\n';
		}
		std::cout << "EOD\n";
		std::cout << "plot sample ";
		for (size_t i = 0; i < spline.Size(); ++i) {
			std::cout << '[' << spline.x[i] << ':' << spline.x[i+1] << "] " << spline[i] << " ls 1, \\\n";
		}
		std::cout << "$points with points ls 1";
	}
}