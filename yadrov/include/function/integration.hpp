#pragma once

#include <functional>
#include <cmath>

#include "polynomial.hpp"
#include "interpolation_polynomial.hpp"

template <class T>
T RectangleMethod(const std::function<T(T)>& f, T a, T b, T h) {
	T res = 0;
	int count = std::floor((b - a) / h);
	T x = a;
	for (int i = 0; i < count; ++i, x += h) {
		res += f(x) * h;
	}
	return res;
}

template <class T>
T TrapeziumMethod(const std::function<T(T)>& f, T a, T b, T h) {
	T res = 0;
	int count = std::floor((b - a) / h);
	T prev_x = a, x = a + h;
	for (int i = 0; i < count; ++i, prev_x = x, x += h) {
		res += (f(prev_x) + f(x)) / 2 * h;
	}
	return res;
}

template <class T>
T SimpsonMethod(const std::function<T(T)>& f, T a, T b, T h) {
	T res = 0;
	int count = std::floor((b - a) / h) / 2;
	T x = a, next_x;
	for (int i = 0; i < count; ++i, x = next_x) {
		next_x = x + h + h;
		Polynomial<T> p = LagrangePolynomial<T>({x, x+h, next_x}, {f(x), f(x+h), f(next_x)});
		res += p.Integrate(x, next_x);
	}
	return res;
}

template <class T>
T RungeRombergError(const std::function<T(std::function<T(T)>, T, T, T)>& IntegrationMethod, const std::function<T(T)>& f, T a, T b, T h, T r, int p) {
	T Ih = IntegrationMethod(f, a, b, h);
	T Irh = IntegrationMethod(f, a, b, r * h);
	return (Ih - Irh) / (std::pow(r, p) - 1);
}

template <class T>
T RungeRombergError(const std::function<T(std::function<T(T)>, T, T, T)>& IntegrationMethod, const std::function<T(T)>& f, T Ih, T a, T b, T h, T r, int p) {
	T Irh = IntegrationMethod(f, a, b, r * h);
	return (Ih - Irh) / (std::pow(r, p) - 1);
}