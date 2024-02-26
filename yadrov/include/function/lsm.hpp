#pragma once

#include <functional>
#include <vector>

#include "../linear/matrix.hpp"
#include "../linear/vector.hpp"
#include "../linear/lup.hpp"

template <class T>
std::vector<T> LSM(const std::vector<T>& x, const std::vector<T>& y, const std::vector<std::function<T(T)>>& basis) {
	if (x.size() != y.size()) {
		throw std::runtime_error("Incompatible arrays");
	}

	Matrix<T> Phi(x.size(), basis.size());
	for (size_t i = 0; i < x.size(); ++i) {
		for (size_t j = 0; j < basis.size(); ++j) {
			Phi[i][j] = basis[j](x[i]);
		}
	}
	Matrix<T> PhiT = Phi.Transpose();
	Matrix<T> G = PhiT * Phi;
	Vector<T> z = PhiT * Vector<T>(y);

	return LUP(G).Solve(z);
}

template <class T>
T SquareError(const std::function<T(T)>& f, const std::vector<T>& x, const std::vector<T>& y) {
	T res = 0;
	for (size_t i = 0; i < x.size(); ++i) {
		res += (y[i] - f(x[i])) * (y[i] - f(x[i]));
	}
	return res;
}