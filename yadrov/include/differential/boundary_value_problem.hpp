#pragma once

#include <functional>
#include <vector>

#include "cauchy_problem.hpp"
#include "../linear/tridiagonal_matrix.hpp"

template <class T>
std::vector<std::vector<T>> ShootingMethod(const function_t<T>& f, const function_t<T>& g, T y1, T y2, T start, T end, T h) {
	std::function<T(T)> Phi = [y2](T y){ return y - y2; };

	T mu0 = 1, mu1 = 2, phi0, phi1, dphi;
	std::vector<std::vector<T>> solution = RungeKuttaMethod(f, g, y1, mu0, start, end, h);
	phi0 = Phi(solution[1].back());
	solution = RungeKuttaMethod(f, g, y1, mu1, start, end, h);
	phi1 = Phi(solution[1].back());

	while (std::abs(phi1) > EPS) {
		dphi = (phi1 - phi0) / (mu1 - mu0);

		mu0 = mu1;
		mu1 -= phi1 / dphi;

		phi0 = phi1;
		solution = RungeKuttaMethod(f, g, y1, mu1, start, end, h);
		phi1 = Phi(solution[1].back());
	}

	return solution;
}

template <class T>
std::vector<std::vector<T>> FiniteDifferenceMethod(const std::function<T(T)>& p, const std::function<T(T)>& q, const std::function<T(T)>& r, T alpha1, T beta1, T gamma1, T alpha2, T beta2, T gamma2, T start, T end, T h) {
	std::vector<std::vector<T>> res(2, std::vector<T>());
	T x = start;
	res[0].push_back(x);

	std::vector<T> a, b, c, d;

	a.push_back(0);
	b.push_back(-alpha1 / h + beta1);
	c.push_back(alpha1 / h);
	d.push_back(gamma1);
	while (x + 2 * h < end || std::abs(end - x - 2 * h) < EPS) {
		x += h; 
		res[0].push_back(x);
		a.push_back(1 / (h * h) - p(x) / (2 * h));
		b.push_back(- 2 / (h * h) + q(x));
		c.push_back(1 / (h * h) + p(x) / (2 * h));
		d.push_back(-r(x));
	}
	x += h; 
	res[0].push_back(x);
	a.push_back(-alpha2 / h);
	b.push_back(alpha2 / h + beta2);
	c.push_back(0);
	d.push_back(gamma2);

	TDMatrix<T> matrix(a, b, c);
	res[1] = matrix.Solve(d);
	return res;
}