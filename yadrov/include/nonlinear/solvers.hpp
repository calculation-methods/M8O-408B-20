#pragma once

#include <functional>
#include <cmath>

template <class T>
using interval_t = std::pair<T, T>;

template <class T>
int sign(T x) {
	if (x > 0) {
		return 1;
	} else if (x < 0) {
		return -1;
	} else {
		return 0;
	}
}

template <class T>
std::pair<T, size_t> IterationMethod(T x, const interval_t<T>& interval, const std::function<T(T)>& f, const std::function<T(T)>& df, double eps) {
	T lambda = sign(df(x)) / std::max(std::abs(df(interval.first)), std::abs(df(interval.second)));
	double q = std::max(std::abs(1 - lambda * df(interval.first)), std::abs(1 - lambda * df(interval.second)));
	q = q / (1 - q);

	double eps_k;
	T next_x;
	size_t iter_count = 0;
	do {
		next_x = x - lambda * f(x);
		eps_k = q * std::abs(next_x - x);
		
		std::swap(next_x, x);
		++iter_count;
	} while (eps_k > eps);

	return {x, iter_count};
}

template <class T>
std::pair<T, size_t> NewtonMethod(T l, T r, const std::function<T(T)>& f, const std::function<T(T)>& df, const std::function<T(T)>& ddf, double eps) {
	T x = l;
	if (f(x) * ddf(x) < eps) {
		x = r;
	}
	double eps_k;
	T next_x;
	size_t iter_count = 0;
	do {
		next_x = x - f(x) / df(x);
		eps_k = std::abs(next_x - x);
		
		std::swap(next_x, x);
		++iter_count;
	} while (eps_k > eps);

	return {x, iter_count};
}