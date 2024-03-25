#pragma once

#include <cmath>
#include <string>

#include "../include/partial_differential/parabolic_pde.hpp"

using std::exp, std::sin, std::cos;
using namespace ParabolicPDE;

const double PI = 2 * std::acos(0.0);

struct Preset {
	std::string name;
	PDE<double> pde;
	std::string f;
	std::string psi;
	std::string gamma1;
	std::string gamma2;
	std::string solution;
	int h_count;
	double t_end;

	Preset() = default;

	Preset(std::string name, PDE<double>&& pde, std::string f, std::string psi, std::string gamma1, std::string gamma2, std::string solution, int h_count, double t_end) :
	       name(name), pde(pde), f(f), psi(psi), gamma1(gamma1), gamma2(gamma2), solution(solution), h_count(h_count), t_end(t_end) {}
};

void GeneratePresets(std::vector<Preset>& presets) {
	presets.resize(3);
	presets[0] = Preset{
		"Sinus",
		PDE<double>{
			1., 0., 0., [](double, double){ return 0; },
			[](double x){ return sin(PI * x); },
			0., 1.,
			0., 1., [](double){ return 0; },
			0., 1., [](double){ return 0; }
		},
		"0", "sin(pi * x)", "0", "0",
		"exp(a pi^2 t) * sin(pi * x)",
		10, 0.25
	};
	presets[0].pde.SetSolution(
		[&pde = presets[0].pde](double x, double t){ return exp(- pde.a * PI * PI * t) * sin(PI * x); });

	presets[1] = Preset{
		"Task 10",
		PDE<double>{
			1., 1., -1., [](double, double){ return 0; }
		},
		"0", "sin(x)", "exp((c-a)t)*(cos(bt)+sin(bt))", "-exp((c-a)t)*(cos(bt)+sin(bt))",
		"exp((c-a)t) * sin(x + bt)",
		20, 1.5
	};
	presets[1].pde.SetBoundaries(
		[](double x){ return sin(x); },
		0., PI,
		1., 1., [&pde = presets[1].pde](double t){ return exp((pde.c - pde.a) * t) * (cos(pde.b * t) + sin(pde.b * t)); },
		1., 1., [&pde = presets[1].pde](double t){ return -exp((pde.c - pde.a) * t) * (cos(pde.b * t) + sin(pde.b * t)); });
	presets[1].pde.SetSolution(
		[&pde = presets[1].pde](double x, double t){ return exp((pde.c - pde.a) * t) * sin(x + pde.b * t); });

	presets[2] = Preset{
		"Task 9",
		PDE<double>{
			1., 1., 0., [](double, double){ return 0; }
		},
		"0", "cos(x)", "-exp(-at)*(cos(bt)+sin(bt))", "exp(-at)*(cos(bt)+sin(bt))",
		"exp((c-a)t) * cos(x + bt)",
		20, 1.5
	};
	presets[2].pde.SetBoundaries(
		[](double x){ return cos(x); },
		0., PI,
		1., -1., [&pde = presets[2].pde](double t){ return -exp(-pde.a * t) * (cos(pde.b * t) + sin(pde.b * t)); },
		1., -1., [&pde = presets[2].pde](double t){ return exp(-pde.a * t) * (cos(pde.b * t) + sin(pde.b * t)); });
	presets[2].pde.SetSolution(
		[&pde = presets[2].pde](double x, double t){ return exp(-pde.a * t) * cos(x + pde.b * t); });
}