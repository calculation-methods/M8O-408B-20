#pragma once

#include <functional>
#include <vector>
#include <stdexcept>

#include "../linear/matrix.hpp"
#include "../linear/vector.hpp"
#include "../linear/lup.hpp"
#include "../function/optimization.hpp"

template <class T, class... Args>
class Matrix<std::function<T(Args...)>> {
private:
	using function_t = std::function<T(Args...)>;

	std::vector<std::vector<function_t>> _data;
	size_t _size;

	Matrix() : _data(), _size(0) { }

public:
	Matrix(size_t n) : _data(n, std::vector<function_t>(n)), _size(n) { }

	Matrix(std::initializer_list<std::vector<function_t>> list) : _data(list), _size(_data.size()) {
		for (const std::vector<function_t>& row: _data) {
			if (row.size() != _size) {
				throw std::runtime_error("Incorrect initializer list");
			}
		}
	}

	std::vector<function_t>& operator[](size_t index) {
		return _data[index];
	}

	const std::vector<function_t>& operator[](size_t index) const {
		return _data[index];
	}

	size_t Size() const {
		return _size;
	}

	Matrix<T> operator()(Args... args) const {
		Matrix<T> res(_size);
		for (size_t i = 0; i < _size; ++i) {
			for (size_t j = 0; j < _size; ++j) {
				res[i][j] = _data[i][j](args...);
			}
		}
		return res;
	}
};

template <class T, class... Args>
class Vector<std::function<T(Args...)>> {
private:
	using function_t = std::function<T(Args...)>;

	std::vector<function_t> _data;
	size_t _size;

	Vector() : _data(), _size(0) { }

public:
	Vector(size_t n) : _data(n), _size(n) { }

	Vector(std::initializer_list<function_t> list) : _data(list), _size(_data.size()) { }

	function_t& operator[](size_t index) {
		return _data[index];
	}

	const function_t& operator[](size_t index) const {
		return _data[index];
	}

	size_t Size() const {
		return _size;
	}

	Vector<T> operator()(Args... args) const {
		Vector<T> res(_size);
		for (size_t i = 0; i < _size; ++i) {
			res[i] = _data[i](args...);
		}
		return res;
	}
};

template <class T>
using function_t = std::function<T(Vector<T>)>;

template <class T>
size_t IterationMethod(Vector<T>& x, const Vector<T>& a, const Vector<T>& b, const Vector<function_t<T>>& F, const Matrix<function_t<T>>& J, double eps) {
	double eps_k;
	Vector<T> next_x(x.Size());
	size_t iter_count = 0;

	Matrix<T> lambda = LUP(J(x)).Invert(), E = Matrix<T>::Identity(x.Size());
	function_t<T> Jphi = [&E, &lambda, &J](Vector<T> x){ return (E - lambda * J(x)).Norm(); };

	std::function<bool(Vector<T>)> region = [&a, &b](Vector<T> x) {
		bool in = true;
		for (size_t i = 0; i < x.Size(); ++i) {
			if (x[i] < a[i] || x[i] > b[i]) {
				in = false;
			}
		}
		return in;
	};
	std::vector<Vector<T>> simplex{a};
	Vector<T> point = a;
	for (size_t i = 0; i < a.Size(); ++i) {
		point[i] = b[i];
		simplex.push_back(point);
		point[i] = a[i];
	}

	T q = Jphi(NelderMeadMethod<T, std::greater<T>>(Jphi, simplex, eps, region));
	if (q > 1 - eps) {
		throw std::runtime_error("Incorrect interval");
	}
	q = q / (1 - q);

	do {
		next_x = x - lambda * F(x);
		eps_k = q * (next_x - x).Norm();

		std::swap(next_x, x);
		++iter_count;
	} while(eps_k > eps);

	return iter_count;
}

template <class T>
size_t NewtonMethod(Vector<T>& x, const Vector<function_t<T>>& F, const Matrix<function_t<T>>& J, double eps, int k) {
	double eps_k;
	Vector<T> dx(x.Size());
	size_t iter_count = 0;
	LUP<T> lu(F.Size());
	if (k == 0) {
		lu.Assign(J(x));
	}

	do {
		if (k > 0 && iter_count % k == 0) {
			lu.Assign(J(x));
		}

		dx = lu.Solve(-F(x));
		x += dx;
		eps_k = dx.Norm();

		++iter_count;
	} while (eps_k > eps);

	return iter_count;
}