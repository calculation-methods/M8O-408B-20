#pragma once

#include <vector>
#include <stdexcept>

template <class T>
T TableFirstDerivative(const std::vector<T>& x, const std::vector<T>& y, T value) {
	if (x.size() != y.size()) {
		throw std::runtime_error("Incompatible arrays");
	}

	if (x.size() < 3) {
		throw std::runtime_error("Too few points");
	}

	if (value < x[0] || value > x.back()) {
		throw std::runtime_error("Out of range");
	}

	size_t i = 0;
	T min = value - x[0];
	for (size_t j = 1; j < x.size(); ++j) {
		if (min > std::abs(value - x[j])) {
			min = std::abs(value - x[j]);
			i = j - 1;
		}
	}

	if (i == x.size() - 2) {
		--i;
	}

	return (y[i+1] - y[i]) / (x[i+1] - x[i]) + ((y[i+2] - y[i+1]) / (x[i+2] - x[i+1]) - (y[i+1] - y[i]) / (x[i+1] - x[i])) / (x[i+2] - x[i]) * (2 * value - x[i] - x[i+1]);
}

template <class T>
T TableSecondDerivative(const std::vector<T>& x, const std::vector<T>& y, T value) {
	if (x.size() != y.size()) {
		throw std::runtime_error("Incompatible arrays");
	}

	if (x.size() < 3) {
		throw std::runtime_error("Too few points");
	}

	if (value < x[0] || value > x.back()) {
		throw std::runtime_error("Out of range");
	}

	size_t i = 0;
	T min = value - x[0];
	for (size_t j = 1; j < x.size(); ++j) {
		if (min > std::abs(value - x[j])) {
			min = std::abs(value - x[j]);
			i = j - 1;
		}
	}

	if (i == x.size() - 2) {
		--i;
	}

	return 2 * ((y[i+2] - y[i+1]) / (x[i+2] - x[i+1]) - (y[i+1] - y[i]) / (x[i+1] - x[i])) / (x[i+2] - x[i]);
}