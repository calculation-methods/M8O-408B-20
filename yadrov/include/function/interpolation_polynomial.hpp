#pragma once

#include <vector>
#include <functional>

#include "polynomial.hpp"

template <class T>
Polynomial<T> LagrangePolynomial(const std::vector<T>& x, const std::vector<T>& y) {
	if (x.size() != y.size()) {
		throw std::runtime_error("Incompatible arrays");
	}
	Polynomial<T> p(x.size() - 1);

	T div;
	Polynomial<T> l;
	for (size_t i = 0; i < x.size(); ++i) {
		l.Assign(0, 1);
		div = 1;

		for (size_t j = 0; j < x.size(); ++j) {
			if (j == i) {
				continue;
			}
			div *= x[i] - x[j];
			l *= Polynomial<T>{-x[j], 1}; // x - x_j
		}
		l *= y[i] / div;
		p += l;
	}

	return p;
}

template <class T>
T DividedDifference(const std::vector<T>& x, const std::vector<T>& y) {
	if (x.size() != y.size()) {
		throw std::runtime_error("Incompatible arrays");
	}

	T res = 0, res_i;
	for (size_t i = 0; i < x.size(); ++i) {
		res_i = 1;
		for (size_t j = 0; j < x.size(); ++j) {
			if (i == j) {
				continue;
			}
			res_i *= x[i] - x[j];
		}

		res += y[i] / res_i;
	}

	return res;
}

template <class T>
Polynomial<T> NewtonPolynomial(const std::vector<T>& x, const std::vector<T>& y) {
	if (x.size() != y.size()) {
		throw std::runtime_error("Incompatible arrays");
	}

	Polynomial<T> p, p_i;
	std::vector<T> x_i, y_i;
	for (size_t i = 0; i < x.size(); ++i) {
		x_i.push_back(x[i]);
		y_i.push_back(y[i]);
		p_i.Assign(0, DividedDifference(x_i, y_i));

		for (size_t j = 0; j < i; ++j) {
			p_i *= Polynomial<T>{-x[j], 1};
		}

		p += p_i;
	}

	return p;
}