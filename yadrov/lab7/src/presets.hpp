#pragma once

#include <map>
#include <cmath>

#include "../../include/partial_differential/elliptic_pde.hpp"

using std::sin, std::cos, std::exp;
using namespace EllipticPDE;

const double pi = 2 * std::acos(0), e = std::exp(1);

void GeneratePresets(std::map<int, PDE<double>>& presets) {
	presets.clear();

	// task 1
	// условия первого рода
	presets[1] = {
		1, 0, 0, 0, [](double, double){ return 0; }, 
		0, 1, 0, 1,
		0, 1, [](double y){ return y; }, 0, 1, [](double y){ return y + 1; },
		0, 1, [](double x){ return x; }, 0, 1, [](double x){ return x + 1; },
		[](double x, double y){ return x + y; }
	};

	// task 8
	// условия первого рода, но больше слагаемых в уравнении
	presets[8] = {
		1, -2, 0, -3, [](double, double){ return 0; }, 
		0, pi/2, 0, pi/2,
		0, 1, [](double y){ return cos(y); }, 0, 1, [](double){ return 0; },
		0, 1, [](double x){ return exp(-x) * cos(x); }, 0, 1, [](double){ return 0; },
		[](double x, double y){ return exp(-x) * cos(x) * cos(y); }
	};

	// task 3
	// условия второго рода с одной стороны
	presets[3] = {
		1, 0, 0, 0, [](double, double){ return 0; }, 
		0, 1, 0, pi/2,
		0, 1, [](double y){ return cos(y); }, 0, 1, [](double y){ return e * cos(y); },
		1, 0, [](double){ return 0; }, 1, 0, [](double x){ return -exp(x); },
		[](double x, double y){ return exp(x) * cos(y); }
	};

	// task 4
	// условия второго рода с другой стороны
	presets[4] = {
		1, 0, 0, 0, [](double, double){ return 0; }, 
		0, pi, 0, 1,
		1, 0, [](double y){ return exp(y); }, 1, 0, [](double y){ return -exp(y); },
		0, 1, [](double x){ return sin(x); }, 0, 1, [](double x){ return e * sin(x); },
		[](double x, double y){ return exp(y) * sin(x); }
	};

	// task 5
	// условия третьего рода
	presets[5] = {
		1, 0, 0, -1, [](double, double){ return 0; },
		0, 1, 0, pi/2,
		1, 0, [](double y){ return cos(y); }, 1, -1, [](double){ return 0; },
		0, 1, [](double x){ return x; }, 0, 1, [](double){ return 0; },
		[](double x, double y){ return x * cos(y); }
	};

	// task 2
	// условия второго рода с разных сторон
	presets[2] = {
		1, 0, 0, 0, [](double, double){ return 0; },
		0, 1, 0, 1,
		1, 0, [](double){ return 0; }, 0, 1, [](double y){ return 1 - y * y; },
		1, 0, [](double){ return 0; }, 0, 1, [](double x){ return x * x - 1; },
		[](double x, double y){ return x * x - y * y; }
	};
}