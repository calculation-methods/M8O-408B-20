#pragma once

#include <functional>
#include <vector>

template <class T>
using function_t = std::function<T(T, T, T)>;

const double EPS = 1e-6;

template <class T>
std::vector<std::vector<T>> EulerMethod(const function_t<T>& f, const function_t<T>& g, T y0, T z0, T start, T end, T h) {
	std::vector<std::vector<T>> res(3, std::vector<T>());
	T x = start, y = y0, z = z0;
	res[0].push_back(x);
	res[1].push_back(y);
	res[2].push_back(z);

	T pred_y, pred_z, next_x, next_y, next_z;
	while (x + h < end || std::abs(end - x - h) < EPS) {
		pred_y = y + h * f(x, y, z);
		pred_z = z + h * g(x, y, z);

		next_x = x + h;
		next_y = y + h * (f(x, y, z) + f(next_x, pred_y, pred_z)) / 2;
		next_z = z + h * (g(x, y, z) + g(next_x, pred_y, pred_z)) / 2;

		x = next_x;
		y = next_y;
		z = next_z;

		res[0].push_back(x);
		res[1].push_back(y);
		res[2].push_back(z);
	}

	return res;
}

template <class T>
std::vector<std::vector<T>> RungeKuttaMethod(const function_t<T>& f, const function_t<T>& g, T y0, T z0, T start, T end, T h) {
	std::vector<std::vector<T>> res(3, std::vector<T>());
	T x = start, y = y0, z = z0;
	res[0].push_back(x);
	res[1].push_back(y);
	res[2].push_back(z);

	T K1, K2, K3, K4, L1, L2, L3, L4;
	while (x + h < end || std::abs(end - x - h) < EPS) {
		K1 = h * f(x, y, z);
		L1 = h * g(x, y, z);

		K2 = h * f(x + h / 2, y + K1 / 2, z + L1 / 2);
		L2 = h * g(x + h / 2, y + K1 / 2, z + L1 / 2);

		K3 = h * f(x + h / 2, y + K2 / 2, z + L2 / 2);
		L3 = h * g(x + h / 2, y + K2 / 2, z + L2 / 2);

		K4 = h * f(x + h, y + K3, z + L3);
		L4 = h * g(x + h, y + K3, z + L3);

		x += h;
		y += (K1 + 2 * K2 + 2 * K3 + K4) / 6;
		z += (L1 + 2 * L2 + 2 * L3 + L4) / 6;

		res[0].push_back(x);
		res[1].push_back(y);
		res[2].push_back(z);
	}

	return res;
}

template <class T>
std::vector<std::vector<T>> AdamsMethod(const function_t<T>& f, const function_t<T>& g, T y0, T z0, T start, T end, T h) {
	std::vector<std::vector<T>> res = RungeKuttaMethod(f, g, y0, z0, start, start + 3 * h, h);
	T x = res[0].back(), y = res[1].back(), z = res[2].back();

	size_t k = 3;
	T pred_y, pred_z, dy, dz;
	while (x + h < end || std::abs(end - x - h) < EPS) {
		pred_y = y + h * (55 * f(res[0][k],   res[1][k],   res[2][k]) - 
		                  59 * f(res[0][k-1], res[1][k-1], res[2][k-1]) + 
		                  37 * f(res[0][k-2], res[1][k-2], res[2][k-2]) - 
		                  9  * f(res[0][k-3], res[1][k-3], res[2][k-3])) / 24;

		pred_z = z + h * (55 * g(res[0][k],   res[1][k],   res[2][k]) - 
		                  59 * g(res[0][k-1], res[1][k-1], res[2][k-1]) + 
		                  37 * g(res[0][k-2], res[1][k-2], res[2][k-2]) - 
		                  9  * g(res[0][k-3], res[1][k-3], res[2][k-3])) / 24;

		dy = h * (9  * f(x + h,       pred_y,      pred_z) + 
		          19 * f(res[0][k],   res[1][k],   res[2][k]) - 
		          5  * f(res[0][k-1], res[1][k-1], res[2][k-1]) + 
		               f(res[0][k-2], res[1][k-2], res[2][k-2])) / 24;

		dz = h * (9  * g(x + h,       pred_y,      pred_z) + 
		          19 * g(res[0][k],   res[1][k],   res[2][k]) - 
		          5  * g(res[0][k-1], res[1][k-1], res[2][k-1]) + 
		               g(res[0][k-2], res[1][k-2], res[2][k-2])) / 24;

		x += h;
		y += dy;
		z += dz;
		k += 1;

		res[0].push_back(x);
		res[1].push_back(y);
		res[2].push_back(z);
	}

	return res;
}

template <class T>
T RungeRombergError(const std::vector<T>& uh, const std::vector<T>& urh, int p) {
	double r = 2;
	double coeff = 1. / (std::pow(r, p) - 1);
	T err = 0;
	for (size_t i = 0; i < urh.size(); ++i) {
		err = std::max(err, std::abs(uh[2 * i] - urh[i]) * coeff);
	}
	return err;
}