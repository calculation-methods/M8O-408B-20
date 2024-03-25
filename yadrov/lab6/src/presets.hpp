#pragma once

#include <cmath>
#include <string>

#include "../include/partial_differential/hyperbolic_pde.hpp"

using std::exp, std::sin, std::cos, std::sqrt;
using namespace HyperbolicPDE;

const double PI = 2 * std::acos(0.0);

struct Preset {
	std::string name;
	PDE<double> pde;
	std::string f;
	std::string psi1;
	std::string psi2;
	std::string gamma1;
	std::string gamma2;
	std::string solution;
	int h_count;
	double t_end;

	Preset() = default;

	Preset(std::string name, PDE<double>&& pde, std::string f, std::string psi1, std::string psi2, std::string gamma1, std::string gamma2, std::string solution, int h_count, double t_end) :
	       name(name), pde(pde), f(f), psi1(psi1), psi2(psi2), gamma1(gamma1), gamma2(gamma2), solution(solution), h_count(h_count), t_end(t_end) {}
};

void GeneratePresets(std::vector<Preset>& presets) {
	presets.resize(5);

	presets[0] = Preset{
		"Sinus",
		PDE<double>{
			1, 0, 0, 0, [](double, double){return 0.;},
			[](double x){return sin(PI * x);}, [](double x){return PI * cos(PI * x);}, [](double x){return -PI * PI * sin(PI * x);}, 
			[](double){return 0.;}, 
			0, 1, 
			0, 1, [](double){return 0.;}, 
			0, 1, [](double){return 0.;}
		},
		"0", "sin(pi x)", "0", "0", "0",
		"cos(pi a t) * sin(pi x)",
		20, 2
	};
	presets[0].pde.SetSolution([&pde = presets[0].pde](double x, double t){return cos(PI * std::sqrt(pde.a) * t) * sin(PI * x);});

	presets[1] = Preset{
		"Task 4",
		PDE<double>{
			1, 0, -5, 0, [](double, double){return 0.;},
			[](double x){return exp(2 * x);}, [](double x){return 2 * exp(2 * x);}, [](double x){return 4 * exp(2 * x);}, 
			[](double){return 0.;}, 
			0, 1, 
			1, -2, [](double){return 0.;}, 
			1, -2, [](double){return 0.;},
			[](double x, double t){return exp(2 * x) * cos(t);}
		},
		"0", "exp(2 x)", "0", "0", "0",
		"exp(2 x) * cos(t)",
		50, PI,
	};

	presets[2] = Preset{
		"Task 2",
		PDE<double>{
			1, 0, 0, 0, [](double, double){return 0.;},
		},
		"0", "sin(x) + cos(x)", "-a (sin(x) + cos(x))", "0", "0",
		"sin(x - at) + cos(x + at)",
		50, 1
	};
	presets[2].pde.SetBoundaries(
		[](double x){return sin(x)+cos(x);}, [](double x){return cos(x)-sin(x);}, [](double x){return -sin(x)-cos(x);},
		[&pde=presets[2].pde](double x){return -pde.a*(sin(x)+cos(x));},
		0, PI,
		1, -1, [](double){return 0.;},
		1, -1, [](double){return 0.;}
	);
	presets[2].pde.SetSolution([&pde=presets[2].pde](double x, double t){return sin(x - sqrt(pde.a) * t) + cos(x + sqrt(pde.a) * t);});
	

	presets[3] = Preset{
		"Task 10",
		PDE<double>{
			1, 1, -1, -3, [](double x, double t){return -cos(x)*exp(-t);},
			[](double x){return sin(x);}, [](double x){return cos(x);}, [](double x){return -sin(x);},
			[](double x){return -sin(x);},
			0, PI,
			1, 0, [](double t){return exp(-t);},
			1, 0, [](double t){return -exp(-t);},
			[](double x, double t){return exp(-t) * sin(x);}
		},
		"-cos(x) * exp(-t)", "sin(x)", "-sin(x)",
		"exp(-t)", "-exp(-t)",
		"exp(-t) * sin(x)",
		100, 2.5
	};

	presets[4] = Preset{
		"Task 8",
		PDE<double>{
			1, 2, -3, -2, [](double, double){return 0;},
			[](double){return 0;}, [](double){return 0;}, [](double){return 0;},
			[](double x){return 2 * exp(-x) * sin(x);},
			0, PI,
			0, 1, [](double){return 0;},
			0, 1, [](double){return 0;},
			[](double x, double t){return exp(-t - x) * sin(x) * sin(2 * t);}
		},
		"-cos(x) * exp(-t)", "sin(x)", "-sin(x)",
		"exp(-t)", "-exp(-t)",
		"exp(-t) * sin(x)",
		50, 3.5
	};
}