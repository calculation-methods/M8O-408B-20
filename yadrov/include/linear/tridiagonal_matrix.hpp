#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "vector.hpp"

template <class T>
class TDMatrix {
private:
	size_t _size;

public:
	std::vector<T> a, b, c;
	
	TDMatrix(size_t n) : _size(n), a(n), b(n), c(n) { }

	TDMatrix(std::vector<T> a, std::vector<T> b, std::vector<T> c) : _size(a.size()), a(a), b(b), c(c) {
		if (a.size() != b.size() || b.size() != c.size()) {
			throw std::runtime_error("Incompatible arrays");
		}
	}

	size_t Size() const {
		return _size;
	}

	Vector<T> Solve(const Vector<T>& vec) {
		if (vec.Size() != _size) {
			throw std::runtime_error("Dimension mismatch");
		}

		Vector<T> p(_size-1), q(_size-1), x(_size);
		p[0] = - c[0] / b[0];
		q[0] = vec[0] / b[0];

		for (size_t i = 1; i < _size-1; ++i) {
			p[i] = - c[i] / (b[i] + a[i] * p[i-1]);
			q[i] = (vec[i] - a[i] * q[i-1]) / (b[i] + a[i] * p[i-1]);
		}

		x[_size-1] = (vec[_size-1] - a[_size-1] * q[_size-2]) / (b[_size-1] + a[_size-1] * p[_size-2]);
		for (int i = _size - 2; i >= 0; --i) {
			x[i] = p[i] * x[i+1] + q[i];
		}

		return x;
	}

	template <class U>
	friend std::ostream& operator<<(std::ostream&, const TDMatrix<U>&);
	template <class U>
	friend std::istream& operator>>(std::istream&, TDMatrix<U>&);
};

template <class T>
std::ostream& operator<<(std::ostream& os, const TDMatrix<T>& matrix) {
	os << matrix.b[0] << ' ' << matrix.c[0] << '\n';
	for (size_t i = 1; i < matrix.Size()-1; ++i) {
		os << std::setw(8) << matrix.a[i] << ' ';
		os << std::setw(8) << matrix.b[i] << ' ';
		os << std::setw(8) << matrix.c[i] << '\n';
	}
	os << matrix.a.back() << ' ' << matrix.b.back() << '\n';
	return os;
}

template <class T>
std::istream& operator>>(std::istream& is, TDMatrix<T>& matrix) {
	is >> matrix.b[0] >> matrix.c[0];
	for (size_t i = 1; i < matrix.Size()-1; ++i) {
		is >> matrix.a[i] >> matrix.b[i] >> matrix.c[i];
	}
	is >> matrix.a.back() >> matrix.b.back();
	return is;
}
