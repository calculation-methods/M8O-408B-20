#pragma once

#include <tuple>
#include <vector>
#include <functional>

enum class ApproxType : int {
	Linear,
	Quadratic,
	Taylor
};

template <class T, template<class> class PDE>
std::tuple<std::vector<T>, std::vector<T>>
GenerateGrid(const PDE<T>& pde, T t_end, int h_count, double sigma, const std::function<T(int, double, T, T, T)>& CourantCondition) {
	int tau_count = CourantCondition(h_count, sigma, t_end, pde.end - pde.start, pde.a);
	std::vector<T> x(h_count + 1), t(tau_count + 1);
	T h = (pde.end - pde.start) / h_count;
	T tau = t_end / tau_count;

	for (int i = 0; i <= h_count; ++i) {
		x[i] = pde.start + h * i;
	}
	for (int k = 0; k <= tau_count; ++k) {
		t[k] = tau * k;
	}

	return {x, t};
}

template <class T>
struct Boundaries {
	struct Coeffs {
		T alpha, beta;
	};

	Coeffs left, right;
};