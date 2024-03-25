#pragma once

#include <vector>
#include <utility>
#include <stdexcept>

#include "matrix.hpp"
#include "vector.hpp"

template <class T>
struct LUP {
	Matrix<T> L;
	Matrix<T> U;
	std::vector<std::pair<size_t, size_t>> P;

	LUP(size_t n): L(n), U(n), P() { }

	LUP(const Matrix<T>& matrix): L(matrix.Size()), U(matrix.Size()), P() {
		matrix.Decompose(L, U, P);
	}

	void Assign(const Matrix<T>& matrix) {
		if (L.Size() != matrix.Size()) {
			L = Matrix<T>(matrix.Size());
			U = Matrix<T>(matrix.Size());
		}
		matrix.Decompose(L, U, P);
	}

	Matrix<T> Compose() {
		return L * U;
	}

	T Det() {
		T res = 1;
		for (size_t i = 0; i < U.Size(); ++i) {
			res *= U[i][i];
		}
		if (P.size() & 1) {
			res = -res;
		}
		return res;
	}

	Vector<T> Solve(Vector<T> b) {
		if (b.Size() != L.Size()) {
			throw std::runtime_error("Dimension mismatch");
		}
		for (const std::pair<size_t, size_t>& p: P) {
			std::swap(b[p.first], b[p.second]);
		}

		Vector<T> x(L.Size()), z(L.Size());
		for (size_t i = 0; i < b.Size(); ++i) {
			z[i] = b[i];
			for (size_t j = 0; j < i; ++j) {
				z[i] -= L[i][j] * z[j];
			}
		}

		for (int i = b.Size()-1; i >= 0; --i) {
			x[i] = z[i];
			for (int j = b.Size()-1; j > i; --j) {
				x[i] -= U[i][j] * x[j];
			}
			x[i] /= U[i][i];
		}

		return x;
	}

	Matrix<T> Invert() {
		Vector<T> e(L.Size());
		Matrix<T> inv(L.Size());

		for (size_t i = 0; i < L.Size(); ++i) {
			e[i] = 1;
			Vector<T> b = Solve(e);
			for (size_t j = 0; j < L.Size(); ++j) {
				inv[j][i] = b[j];
			}
			e[i] = 0;
		}
		
		return inv;
	}
};
