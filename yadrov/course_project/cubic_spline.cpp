#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>

#include "../include/function/spline.hpp"

/*
-p - ignore p values in input
-plot - gnuplot script
*/

int main(int argc, char* argv[]) {
	bool plot = false;
	bool ignore_p = false;

	for (int i = 1; i < argc; ++i) {
		if (!plot && std::strcmp(argv[i], "-plot") == 0) {
			plot = true;

		} else if (!ignore_p && std::strcmp(argv[i], "-p") == 0) {
			ignore_p = true;

		} else {
			std::cerr << "Unknown option: " << argv[i] << std::endl;
		}
	}

	std::vector<double> x, y;
	double a, b, p;
	while (std::cin >> a >> b) {
		x.push_back(a);
		y.push_back(b);
		if (ignore_p) {
			std::cin >> p;
		}
	}

	CubicSpline<double> spline(x, y);

	if (!plot) {
		std::cout << spline;

	} else {
		std::cout << "set term pdfcairo\nset size ratio -1\nset key off\n";
		const auto [xmin, xmax] = std::minmax_element(spline.x.begin(), spline.x.end());
		const auto [ymin, ymax] = std::minmax_element(spline.y.begin(), spline.y.end());
		double dx = (*xmax - *xmin) * 0.1, dy = (*ymax - *ymin) * 0.5;
		std::cout << "set xrange [" << (*xmin-dx) << ':' << (*xmax+dx) << "]\nset yrange [" << (*ymin-dy) << ':' << (*ymax+dy) << "]\n";
		std::cout << "set style line 1 lw 1.5 lt 4 pt 7 ps 0.5\n";
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