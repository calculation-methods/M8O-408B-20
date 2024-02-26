#pragma once

#include <cmath>
#include <vector>
#include <complex>

#include "matrix.hpp"
#include "vector.hpp"

template <class T>
void QR_Decomposition(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R) {
	R = A;
	Q = Matrix<T>::Identity(A.Size());
	Vector<T> v(A.Size());
	Matrix<T> H(A.Size());
	for (size_t k = 0; k < A.Size() - 1; ++k) {
		for (size_t i = 0; i < k; ++i) {
			v[i] = 0;
		}
		v[k] = R[k][k] * R[k][k];
		for (size_t i = k+1; i < A.Size(); ++i) {
			v[i] = R[i][k];
			v[k] += R[i][k] * R[i][k];
		}
		v[k] = std::sqrt(v[k]) * (R[k][k] >= 0 ? 1 : -1) + R[k][k];

		H = Matrix<T>::Householder(A.Size(), v);
		R = H * R;
		Q = Q * H;
	}
}

template <class T>
size_t QR_Eigenvalues(Matrix<T> A, double eps, std::vector<std::complex<T>>& eigenvalues) {
	Matrix<T> Q(A.Size()), R(A.Size());
	size_t iter_count = 0;
	Vector<double> eps_1(A.Size()), eps_2(A.Size()), eps_3(A.Size());
	eigenvalues.assign(A.Size(), std::complex<T>());
	std::vector<std::complex<T>> prev_eigenvalues(A.Size(), 1e18);

	do {
		QR_Decomposition(A, Q, R);
		A = R * Q;

		// вычисление собственных значений
		eigenvalues.swap(prev_eigenvalues);
		for (size_t i = 0; i < A.Size(); ++i) {
			if (i == A.Size()-1 || std::abs(A[i+1][i]) < eps) {
				eigenvalues[i] = std::complex(A[i][i]);
			} else {
				T d = (A[i][i] + A[i+1][i+1]) * (A[i][i] + A[i+1][i+1]) - 4 * (A[i][i] * A[i+1][i+1] - A[i][i+1] * A[i+1][i]);
				if (d > eps) {
					continue;
				}
				T re = (A[i][i] + A[i+1][i+1]) * 0.5;
				T im = std::sqrt(std::abs(d)) * 0.5;
				if (std::abs(im) < eps) {
					continue;
				}
				eigenvalues[i] = std::complex(re, im);
				eigenvalues[i+1] = std::complex(re, -im);
				++i;
			}
		}
		
		// вычисление погрешностей
		for (size_t j = 0; j < A.Size(); ++j) {
			eps_1[j] = 0;
			for (size_t i = j+1; i < A.Size(); ++i) {
				eps_1[j] += A[i][j] * A[i][j]; 
			}
			eps_1[j] = std::sqrt(eps_1[j]);

			eps_2[j] = 0;
			for (size_t i = j+2; i < A.Size(); ++i) {
				eps_2[j] += A[i][j] * A[i][j]; 
			}
			eps_2[j] = std::sqrt(eps_2[j]);

			eps_3[j] = std::abs(prev_eigenvalues[j]) - std::abs(eigenvalues[j]);
		}

		++iter_count;
	} while ((eps_1.Norm(2) > eps && eps_2.Norm(2) > eps) || eps_3.Norm() > eps);

	return iter_count;
}
