#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <utility>

#include "vector.hpp"

template <class T>
class Matrix {
private:
	std::vector<std::vector<T>> _data;
	std::pair<size_t, size_t> _size;

	Matrix() : _data(), _size{0, 0} { }

public:
	Matrix(size_t n) : _data(n, std::vector<T>(n)), _size{n, n} { }

	Matrix(size_t n, size_t m) : _data(n, std::vector<T>(m)), _size{n, m} { }

	Matrix(std::pair<size_t, size_t> size) : _data(size.first, std::vector<T>(size.second)), _size(size) { }

	Matrix(std::initializer_list<std::vector<T>> list) : _data(list), _size(_data.size(), _data.back().size()) {
		for (const std::vector<T>& row: _data) {
			if (row.size() != _size.second()) {
				throw std::runtime_error("Incorrect initializer list");
			}
		}
	}

	std::vector<T>& operator[](size_t index) {
		return _data[index];
	}

	const std::vector<T>& operator[](size_t index) const {
		return _data[index];
	}

	bool Square() const {
		return _size.first == _size.second;
	}

	size_t Size() const {
		if (_size.first != _size.second) {
			throw std::runtime_error("Non-square matrix");
		}
		return _size.first;
	}

	std::pair<size_t, size_t> Sizes() const {
		return _size;
	}

	size_t Rows() const {
		return _size.first;
	}

	size_t Columns() const {
		return _size.second;
	}

	static Matrix<T> Identity(size_t n) {
		Matrix<T> result(n);
		for (size_t i = 0; i < n; ++i) {
			result[i][i] = 1;
		}
		return result;
	}

	static Matrix<T> Rotation(size_t n, size_t k, size_t l, double angle) {
		Matrix<T> result = Matrix<T>::Identity(n);
		if (k > l) {
			std::swap(k, l);
		}
		double c = std::cos(angle), s = std::sin(angle);
		result[k][k] = c;
		result[l][l] = c;
		result[k][l] = -s;
		result[l][k] = s;
		return result;
	}

	static Matrix<T> Householder(size_t n, const Vector<T>& v) {
		Matrix<T> vv(n);
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				vv[i][j] = v[i] * v[j];
			}
		}

		return Matrix<T>::Identity(n) - vv * (2.0 / (v * v));
	}

	T Norm() const {
		T norm = 0;
		for (const std::vector<T>& row: _data) {
			T sum = 0;
			for (const T& a: row) {
				sum += std::abs(a);
			}
			if (sum > norm) {
				norm = sum;
			}
		}
		return norm;
	}

	Matrix<T> Transpose() {
		Matrix<T> result(std::make_pair(_size.second, _size.first));
		for (size_t i = 0; i < _size.second; ++i) {
			for (size_t j = 0; j < _size.first; ++j) {
				result[i][j] = _data[j][i];
			}
		}
		return result;
	}

	void SwapRows(size_t i, size_t j) {
		if (i >= _size.first || j >= _size.first) {
			throw std::runtime_error("Invalid index");
		}

		std::swap(_data[i], _data[j]);
	}

	void SwapRows(std::vector<std::pair<size_t, size_t>>& P) {
		for (const std::pair<size_t, size_t>& p: P) {
			if (p.first >= _size.first || p.second >= _size.first) {
				throw std::runtime_error("Invalid index");
			}
		}

		for (const std::pair<size_t, size_t>& p: P) {
			std::swap(_data[p.first], _data[p.second]);
		}
	}

	void SwapColumns(size_t i, size_t j) {
		if (i >= _size.second || j >= _size.second) {
			throw std::runtime_error("Invalid index");
		}

		for (size_t k = 0; k < _size.first; ++k) {
			std::swap(_data[k][i], _data[k][j]);
		}
	}

	bool Decompose(Matrix<T>& L, Matrix<T>& U, std::vector<std::pair<size_t, size_t>>& P) const {
		if (Size() != L.Size() || Size() != U.Size()) {
			throw std::runtime_error("Dimension mismatch");
		}
		P.clear();
		for (size_t i = 0; i < Size(); ++i) {
			for (size_t j = 0; j < Size(); ++j) {
				L[i][j] = 0;
				U[i][j] = _data[i][j];
			}
			L[i][i] = 1;
		}

		for (size_t k = 0; k < Size()-1; ++k) {
			T max = U[k][k];
			size_t max_ind = k;
			for (size_t i = k+1; i < Size(); ++i) {
				if (U[i][k] > max) {
					max = U[i][k];
					max_ind = i;
				}
			}

			if (max_ind != k) {
				P.push_back(std::make_pair(k, max_ind));
				U.SwapRows(k, max_ind);
				L.SwapRows(k, max_ind);
				L.SwapColumns(k, max_ind);
			}

			for (size_t i = k+1; i < Size(); ++i) {
				L[i][k] = U[i][k] / U[k][k];

				for (size_t j = 0; j < Size(); ++j) {
					U[i][j] -= L[i][k] * U[k][j];
				}
			}
		}

		return true; 
	}
};

template <class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
	for (size_t i = 0; i < matrix.Rows(); ++i) {
		for (size_t j = 0; j < matrix.Columns(); ++j) {
			os << std::setw(8) << matrix[i][j] << ' ';
		}
		os << '\n';
	}
	return os;
}

template <class T>
std::istream& operator>>(std::istream& is, Matrix<T>& matrix) {
	for (size_t i = 0; i < matrix.Rows(); ++i) {
		for (size_t j = 0; j < matrix.Columns(); ++j) {
			is >> matrix[i][j];
		}
	}
	return is;
}

template <class T, class U>
Matrix<T> operator*(const Matrix<T>& a, const Matrix<U>& b) {
	if (a.Columns() != b.Rows()) {
		throw std::runtime_error("Dimension mismatch");
	}
	Matrix<T> res(std::make_pair(a.Rows(), b.Columns()));
	for (size_t i = 0; i < a.Rows(); ++i) {
		for (size_t j = 0; j < b.Columns(); ++j) {
			for (size_t k = 0; k < a.Columns(); ++k) {
				res[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return res;
}

template <class T, class U>
Matrix<T> operator+(const Matrix<T>& a, const Matrix<U>& b) {
	if (a.Sizes() != b.Sizes()) {
		throw std::runtime_error("Dimension mismatch");
	}
	Matrix<T> res(a.Sizes());
	for (size_t i = 0; i < a.Rows(); ++i) {
		for (size_t j = 0; j < a.Columns(); ++j) {
			res[i][j] = a[i][j] + b[i][j];
		}
	}
	return res;
}

template <class T, class U>
Matrix<T> operator-(const Matrix<T>& a, const Matrix<U>& b) {
	if (a.Sizes() != b.Sizes()) {
		throw std::runtime_error("Dimension mismatch");
	}
	Matrix<T> res(a.Sizes());
	for (size_t i = 0; i < a.Rows(); ++i) {
		for (size_t j = 0; j < a.Columns(); ++j) {
			res[i][j] = a[i][j] - b[i][j];
		}
	}
	return res;
}

template <class T, class U>
Vector<T> operator*(const Matrix<T>& a, const Vector<U>& b) {
	if (a.Columns() != b.Size()) {
		throw std::runtime_error("Dimension mismatch");
	}
	Vector<T> res(a.Rows());
	for (size_t i = 0; i < a.Rows(); ++i) {
		for (size_t j = 0; j < b.Size(); ++j) {
			res[i] += a[i][j] * b[j];
		}
	}
	return res;
}

template <class T, class U>
Matrix<T> operator*(const Matrix<T>& a, const U& b) {
	Matrix<T> res(a.Sizes());
	for (size_t i = 0; i < a.Rows(); ++i) {
		for (size_t j = 0; j < a.Columns(); ++j) {
			res[i][j] = a[i][j] * b;
		}
	}
	return res;
}