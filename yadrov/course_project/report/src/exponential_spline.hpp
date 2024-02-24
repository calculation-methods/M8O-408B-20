#pragma once

#include <vector>
#include <cmath>
#include <iostream>

#include "../linear/tridiagonal_matrix.hpp"
#include "../linear/vector.hpp"

template <class T>
class ExponentialSpline {
public:
	class Segment {
	public:
		T p, h, t1, t2, x1, x2, f1, f2;

		Segment(): p(0), h(0), t1(0), t2(0), x1(0), x2(0), f1(0), f2(0) { }

		Segment(T p, T x1, T x2, T f1, T f2, T t1, T t2):
			p(p), h(std::abs(x1 - x2)), t1(t1), t2(t2), x1(x1), x2(x2), f1(f1), f2(f2) { }

		T operator()(T x) {
			return (t1 * std::sinh(p * (x2 - x)) + t2 * std::sinh(p * (x - x1))) / (p * p * std::sinh(p * h)) + 
			       (f1 - t1 / (p * p)) * (x2 - x) / h +
			       (f2 - t2 / (p * p)) * (x - x1) / h;
		}

		friend std::ostream& operator<<(std::ostream& os, const Segment& seg) {
			os << "(" << seg.t1 << " * sinh(" << seg.p << " * (" << seg.x2 << " - x)) + " << seg.t2 << " * sinh(" << seg.p << " * (x - " << seg.x1 << "))) / " << seg.p * seg.p * std::sinh(seg.p * seg.h) << " + " 
			   << (seg.f1 - seg.t1 / (seg.p * seg.p)) << " * (" << seg.x2 << " - x) / " << seg.h << " + "
			   << (seg.f2 - seg.t2 / (seg.p * seg.p)) << " * (x - " << seg.x1 << ") / " << seg.h;
			return os;
		}
	};

	std::vector<T> x, y, p;
	std::vector<T> h, d, e;
	Vector<T> b, t;
	std::vector<Segment> segment;

	ExponentialSpline(const std::vector<T>& x, const std::vector<T>& y, const std::vector<T>& p) : x(x), y(y), p(p), h(x.size() - 1), d(x.size() - 1), e(x.size() - 1), b(x.size()), t(x.size()), segment(x.size() - 1) {
		if (x.size() != y.size() && x.size() - 1 != p.size()) {
			throw std::runtime_error("Incompatible arrays");
		}

		Solve();
	}

	void Solve() {
		size_t n = h.size();
		for (size_t i = 0; i < x.size() - 1; ++i) {
			h[i] = std::abs(x[i + 1] - x[i]);
			d[i] = (p[i] * std::cosh(p[i] * h[i]) / std::sinh(p[i] * h[i]) - 1 / h[i]) / (p[i] * p[i]);
			e[i] = (1 / h[i] - p[i] / std::sinh(p[i] * h[i])) / (p[i] * p[i]);
		}

		b[0] = 0;
		b[n] = 0;
		for (size_t i = 1; i < n; ++i) {
			b[i] = (y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1];
		}

		TDMatrix<T> matrix(n+1);
		matrix.b[0] = 1;
		for (size_t i = 1; i < n; ++i) {
			matrix.a[i] = e[i-1];
			matrix.b[i] = d[i-1] + d[i];
			matrix.c[i] = e[i];
		}
		matrix.b[n] = 1;

		t = matrix.Solve(b);

		for (size_t i = 0; i < n; ++i) {
			segment[i] = Segment(p[i], x[i], x[i+1], y[i], y[i+1], t[i], t[i+1]);
		}
	}

	void Tense(T relax) {
		for (size_t k = 1; k < t.Size()-1; ++k) {
			if (t[k] * b[k] < 0) {
				T lambda = std::max(std::abs(b[k]), (d[k-1] + d[k]) * std::abs(t[k])) / 2 / std::max(std::abs(t[k-1]), std::abs(t[k+1]));
				for (size_t i = k-1; i <= k; ++i) {
					T new_p = 1 / std::sqrt(lambda * h[i]);
					p[i] = p[i] + relax * (new_p - p[i]);
				}
			}
		}

		Solve();
	}

	size_t Size() const {
		return segment.size();
	}

	const Segment& operator[](size_t index) {
		return segment[index];
	}

	T operator()(const T& value) {
		if (value < x[0] || value > x.back()) {
			throw std::runtime_error("Out of range");
		}

		for (size_t i = 0; i < x.size()-1; ++i) {
			if (value >= x[i] && value <= x[i+1] ) {
				return segment[i](value);
			}
		}

		return 0; // just for silence the warning
	}
};

template <class T>
std::ostream& operator<<(std::ostream& os, const ExponentialSpline<T>& spline) {
	for (size_t i = 0; i < spline.x.size() - 1; ++i) {
		os << '[' << spline.x[i] << ':' << spline.x[i+1] << "] " << spline.segment[i] << '\n';
	}
	return os;
}