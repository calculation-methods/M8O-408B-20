#pragma once

#include <cmath>
#include <stdexcept>

#include "matrix.hpp"
#include "vector.hpp"

template <class T>
size_t IterationMethod(const Matrix<T>& A, const Vector<T>& b, const Matrix<T>& M, Vector<T>& x, double eps) {
	if ((A.Size() != b.Size()) || (b.Size() != M.Size()) || (M.Size() != x.Size())) {
		throw std::runtime_error("Dimension mismatch");
	}
	Matrix<T> alpha = Matrix<T>::Identity(A.Size()) + M * A;
	Vector<T> beta = -(M * b);
	x = beta;

	T alpha_norm = alpha.Norm();
	size_t count = 0;
	double eps_k;
	Vector<T> next_x(x.Size());

	do {
		next_x = alpha * x + beta;
		if (alpha_norm < 1) {
			eps_k = alpha_norm / (1.0 - alpha_norm) * (next_x - x).Norm();
		} else {
			eps_k = (next_x - x).Norm();
		}
		std::swap(next_x, x);
		++count;
	} while (eps_k > eps);

	return count;
}

template <class T>
size_t JacobiMethod(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x, double eps) {
	Matrix<T> M(A.Size());
	for (size_t i = 0; i < A.Size(); ++i) {
		M[i][i] = -1.0 / A[i][i];
	}
	return IterationMethod(A, b, M, x, eps);
}

template <class T>
size_t SeidelMethod(const Matrix<T>& A, const Vector<T>& b, const Matrix<T>& M, Vector<T>& x, double eps) {
	if ((A.Size() != b.Size()) || (b.Size() != M.Size()) || (M.Size() != x.Size())) {
		throw std::runtime_error("Dimension mismatch");
	}
	Matrix<T> alpha = Matrix<T>::Identity(A.Size()) + M * A;
	Vector<T> beta = -(M * b);
	x = beta;

	Matrix<T> alpha2 = alpha;
	for (size_t i = 0; i < alpha2.Size(); ++i) {
		for (size_t j = 0; j < i; ++j) {
			alpha2[i][j] = 0;
		}
	}

	T alpha_norm = alpha.Norm();
	T alpha2_norm = alpha2.Norm();
	size_t count = 0;
	double eps_k;
	Vector<T> next_x(x.Size());

	do {
		for (size_t i = 0; i < x.Size(); ++i) {
			next_x[i] = beta[i];
			for (size_t j = 0; j < i; ++j) {
				next_x[i] += alpha[i][j] * next_x[j];
			}
			for (size_t j = i; j < x.Size(); ++j) {
				next_x[i] += alpha[i][j] * x[j];
			}
		}

		if (alpha_norm < 1) {
			eps_k = alpha2_norm / (1.0 - alpha_norm) * (next_x - x).Norm();
		} else {
			eps_k = (next_x - x).Norm();
		}
		std::swap(next_x, x);
		++count;
	} while (eps_k > eps);

	return count;
}

template <class T>
size_t JacobiSeidelMethod(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x, double eps) {
	Matrix<T> M(A.Size());
	for (size_t i = 0; i < A.Size(); ++i) {
		M[i][i] = -1.0 / A[i][i];
	}
	return SeidelMethod(A, b, M, x, eps);
}
