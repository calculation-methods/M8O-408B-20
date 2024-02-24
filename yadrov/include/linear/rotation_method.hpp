#pragma once

#include <cmath>

#include "matrix.hpp"

template <class T>
size_t RotationMethod(Matrix<T> A, double eps, std::vector<T>& values, Matrix<T>& U) {
	values.assign(A.Size(), 0);
	U = Matrix<T>::Identity(A.Size());

	size_t k, l;
	T max;
	Matrix<T> transform(A.Size());
	double angle, eps_k;
	size_t count = 0;
	do {
		max = 0;
		for (size_t i = 1; i < A.Size(); ++i) {
			for (size_t j = 0; j < i; ++j) {
				if (std::abs(A[i][j]) > max) {
					max = std::abs(A[i][j]);
					l = i;
					k = j;
				}
			}
		}

		angle = 0.5 * std::atan2(2 * A[k][l], A[k][k] - A[l][l]);
		transform = Matrix<T>::Rotation(A.Size(), k, l, angle);

		A = transform.Transpose() * A * transform;
		U = U * transform;

		eps_k = 0;
		for (size_t i = 1; i < A.Size(); ++i) {
			for (size_t j = 0; j < i; ++j) {
				eps_k += A[i][j] * A[i][j];
			}
		}
		eps_k = std::sqrt(eps_k);

		++count;
	} while (eps_k > eps);

	for (size_t i = 0; i < A.Size(); ++i) {
		values[i] = A[i][i];
	}

	return count;
}