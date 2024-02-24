#pragma once

#include <functional>
#include <vector>
#include <tuple>

#include "../linear/tridiagonal_matrix.hpp"
#include "../linear/vector.hpp"
#include "common.hpp"

namespace HyperbolicPDE {
	template <class T>
	using grid_t = std::vector<std::vector<T>>;

	template <class T>
	struct PDE {
		using f_t = std::function<T(T)>;
		using f_x = f_t;
		using f_x_t = std::function<T(T, T)>;

		T a, b, c, d;
		f_x_t f;
		f_x psi1, dpsi1, d2psi1, psi2;
		T start, end;
		T alpha1, beta1;
		f_t gamma1;
		T alpha2, beta2;
		f_t gamma2;
		f_x_t solution;

		PDE() = default;

		PDE(T a, T b, T c, T d, f_x_t f, f_x psi1, f_x dpsi1, f_x d2psi1, f_x psi2, T start, T end, 
		    T alpha1, T beta1, f_t gamma1, T alpha2, T beta2, f_t gamma2, f_x_t solution) : 
		    a(a), b(b), c(c), d(d), f(f), psi1(psi1), dpsi1(dpsi1), d2psi1(d2psi1), psi2(psi2), start(start), end(end), alpha1(alpha1), beta1(beta1), gamma1(gamma1), 
		    alpha2(alpha2), beta2(beta2), gamma2(gamma2), solution(solution) {}

		PDE(T a, T b, T c, T d, f_x_t f, f_x psi1, f_x dpsi1, f_x d2psi1, f_x psi2, T start, T end, 
		    T alpha1, T beta1, f_t gamma1, T alpha2, T beta2, f_t gamma2) : 
		    a(a), b(b), c(c), d(d), f(f), psi1(psi1), dpsi1(dpsi1), d2psi1(d2psi1), psi2(psi2), start(start), end(end), alpha1(alpha1), beta1(beta1), gamma1(gamma1), 
		    alpha2(alpha2), beta2(beta2), gamma2(gamma2) {}

		PDE(T a, T b, T c, T d, f_x_t f) : 
		    a(a), b(b), c(c), d(d), f(f) {}

		void SetEquation(T a_, T b_, T c_, T d_, f_x_t f_) {
			a = a_;
			b = b_;
			c = c_;
			d = d_;
			f = f_;
		}

		void SetBoundaries(f_x psi1_, f_x dpsi1_, f_x d2psi1_, f_x psi2_, T start_, T end_, T alpha1_, T beta1_, std::function<T(T)> gamma1_, 
			            T alpha2_, T beta2_, std::function<T(T)> gamma2_) {
			psi1 = psi1_;
			dpsi1 = dpsi1_;
			d2psi1 = d2psi1_;
			psi2 = psi2_;
			start = start_;
			end = end_;
			alpha1 = alpha1_;
			beta1 = beta1_;
			gamma1 = gamma1_;
			alpha2 = alpha2_;
			beta2 = beta2_;
			gamma2 = gamma2_;
		}

		void SetSolution(f_x_t solution_) {
			solution = solution_;
		}
	};

	template <class T>
	int CourantCondition(int h_count, double sigma, T t_end, T end, T a) {
		return t_end * a * h_count / (end * sigma);
	}

	template <class T>
	void StartConditions(const PDE<T>& pde, const std::vector<T>& x, const std::vector<T>& t, grid_t<T>& u, T tau, ApproxType type) {
		for (size_t i = 0; i < x.size(); ++i) {
			u[0][i] = pde.psi1(x[i]);
			if (type == ApproxType::Linear) {
				u[1][i] = u[0][i] + tau * pde.psi2(x[i]);
			} else { // Taylor
				u[1][i] = tau * (1 + tau * pde.d / 2) * pde.psi2(x[i]) + u[0][i] + tau * tau / 2 * (pde.a * pde.d2psi1(x[i]) + 
				          pde.b * pde.dpsi1(x[i]) + pde.c * pde.psi1(x[i]) + pde.f(x[i], t[0]));
			}
		}
	}

	template <class T>
	Boundaries<T> BoundariesConditions(const PDE<T>& pde, const std::vector<T>& x, const std::vector<T>& t, grid_t<T>& u, T h, T tau, ApproxType type) {
		Boundaries<T> bound;

		if (type == ApproxType::Linear) {
			bound.left.alpha = -pde.alpha1 / h + pde.beta1;
			bound.left.beta = pde.alpha1 / h;

			bound.right.alpha = pde.alpha2 / h + pde.beta2;
			bound.right.beta = -pde.alpha2 / h;

		} else if (type == ApproxType::Quadratic) {
			bound.left.alpha = -3 * pde.alpha1 / (2 * h) + pde.beta1;
			bound.left.beta = 2 * pde.alpha1 / h;

			bound.right.alpha = 3 * pde.alpha2 / (2 * h) + pde.beta2;
			bound.right.beta = -2 * pde.alpha2 / h;

		} else if (type == ApproxType::Taylor) {
			bound.left.alpha = pde.alpha1 * (-1 - h * h / (2 * pde.a) * (1 / (tau * tau) - pde.c - pde.d / tau)) + pde.beta1 * h * (1 - pde.b * h / (2 * pde.a));
			bound.left.beta = pde.alpha1;

			bound.right.alpha = pde.alpha2 * (-1 - h * h / (2 * pde.a) * (1 / (tau * tau) - pde.c - pde.d / tau)) + pde.beta2* h * (-1 - pde.b * h / (2 * pde.a));
			bound.right.beta = pde.alpha2;
		}

		return bound;
	}

	template <class T>
	std::tuple<std::vector<T>, std::vector<T>, grid_t<T>>
	ExplicitSolver(const PDE<T>& pde, T t_end, int h_count, double sigma, ApproxType start_type, ApproxType bound_type) {
		auto [x, t] = GenerateGrid<T, PDE>(pde, t_end, h_count, sigma, CourantCondition<T>);
		return {x, t, ExplicitSolver(pde, x, t, t_end, start_type, bound_type)};
	}

	template <class T>
	grid_t<T> ExplicitSolver(const PDE<T>& pde, const std::vector<T>& x, const std::vector<T>& t, T t_end, ApproxType start_type, ApproxType bound_type) {
		int h_count = x.size() - 1, tau_count = t.size() - 1;
		grid_t<T> u(tau_count + 1, std::vector<T>(h_count + 1));
		T h = (pde.end - pde.start) / h_count;
		T tau = t_end / tau_count;

		StartConditions(pde, x, t, u, tau, start_type);
		Boundaries bound = BoundariesConditions(pde, x, t, u, h, tau, bound_type);

		T gamma1, gamma2, coeff = (1 / tau - pde.d / 2) / tau;
		for (int k = 1; k < tau_count; ++k) {
			for (int i = 1; i < h_count; ++i) {
				T ddu = (u[k][i-1] - 2 * u[k][i] + u[k][i+1]) / (h * h);
				T du = (u[k][i+1] - u[k][i-1]) / (2 * h);
				u[k+1][i] = (pde.a * ddu + pde.b * du + pde.c * u[k][i] + pde.f(x[i], t[k]) - pde.d / (2 * tau) * u[k-1][i] + (2 * u[k][i] - u[k-1][i]) / (tau * tau)) / coeff;
			}

			gamma1 = pde.gamma1(t[k+1]);
			gamma2 = pde.gamma2(t[k+1]);
			if (bound_type == ApproxType::Taylor) {
				gamma1 = gamma1 * h * (1 - pde.b * h / (2 * pde.a)) + pde.alpha1 * h * h / (2 * pde.a) * (-pde.f(x[0], t[k+1]) + pde.d / tau * u[k][0] + (u[k-1][0] - 2 * u[k][0]) / (tau * tau));
				gamma2 = gamma2 * h * (-1 - pde.b * h / (2 * pde.a)) + pde.alpha2 * h * h / (2 * pde.a) * (-pde.f(x[h_count], t[k+1]) + pde.d / tau * u[k][h_count] + (u[k-1][h_count] - 2 * u[k][h_count]) / (tau * tau));
			}

			u[k+1][0] = (gamma1 - u[k+1][1] * bound.left.beta) / bound.left.alpha;
			u[k+1][h_count] = (gamma2 - u[k+1][h_count-1] * bound.right.beta) / bound.right.alpha;

			if (bound_type == ApproxType::Quadratic) {
				u[k+1][0] -= -pde.alpha1 / (2 * h) * u[k+1][2] / bound.left.alpha;
				u[k+1][h_count] -= pde.alpha2 / (2 * h) * u[k+1][h_count-2] / bound.right.alpha;
			}
		}
		return u;
	}

	template <class T>
	std::tuple<std::vector<T>, std::vector<T>, grid_t<T>>
	ImplicitSolver(const PDE<T>& pde, T t_end, int h_count, double sigma, ApproxType start_type, ApproxType bound_type) {
		auto [x, t] = GenerateGrid<T, PDE>(pde, t_end, h_count, sigma, CourantCondition<T>);
		return {x, t, ImplicitSolver(pde, x, t, t_end, start_type, bound_type)};
	}

	template <class T>
	grid_t<T> ImplicitSolver(const PDE<T>& pde, const std::vector<T>& x, const std::vector<T>& t, T t_end, ApproxType start_type, ApproxType bound_type) {
		int h_count = x.size() - 1, tau_count = t.size() - 1;
		grid_t<T> u(tau_count + 1, std::vector<T>(h_count + 1));
		T h = (pde.end - pde.start) / h_count;
		T tau = t_end / tau_count;

		StartConditions(pde, x, t, u, tau, start_type);
		Boundaries bound = BoundariesConditions(pde, x, t, u, h, tau, bound_type);

		T alpha = (pde.a / h - pde.b / 2) / h,
		  beta = - 2 * pde.a / (h * h) + pde.c + (pde.d / 2 - 1 / tau) / tau,
		  gamma = (pde.a / h + pde.b / 2) / h;

		TDMatrix<T> matrix(h_count+1);
		for (int i = 1; i < h_count; ++i) {
			matrix.a[i] = alpha;
			matrix.b[i] = beta;
			matrix.c[i] = gamma;
		}

		Vector<T> v(h_count+1);
		for (int k = 1; k < tau_count; ++k) {
			for (int i = 1; i < h_count; ++i) {
				v[i] = (- 2 * u[k][i] + u[k-1][i]) / (tau * tau) - pde.f(x[i], t[k+1]) + pde.d * u[k-1][i] / (2 * tau);
			}
			v[0] = pde.gamma1(t[k+1]);
			v[h_count] = pde.gamma2(t[k+1]);

			matrix.b[0] = bound.left.alpha;
			matrix.c[0] = bound.left.beta;
			matrix.a[h_count] = bound.right.beta;
			matrix.b[h_count] = bound.right.alpha;
			
			if (bound_type == ApproxType::Quadratic) {
				T coeff = -pde.alpha1 / (2 * h) / gamma;
				matrix.b[0] -= coeff * alpha;
				matrix.c[0] -= coeff * beta;
				v[0] -= coeff * v[1]; 
				
				coeff = pde.alpha2 / (2 * h) / alpha;
				matrix.a[h_count] -= coeff * beta;
				matrix.b[h_count] -= coeff * gamma;
				v[h_count] -= coeff * v[h_count-1];
			
			} else if (bound_type == ApproxType::Taylor) {
				v[0] = v[0] * h * (1 - pde.b * h / (2 * pde.a)) + pde.alpha1 * h * h / (2 * pde.a) * (-pde.f(x[0], t[k+1]) + pde.d / tau * u[k][0] + (u[k-1][0] - 2 * u[k][0]) / (tau * tau));
				v[h_count] = v[h_count] * h * (-1 - pde.b * h / (2 * pde.a)) + pde.alpha2 * h * h / (2 * pde.a) * (-pde.f(x[h_count], t[k+1]) + pde.d / tau * u[k][h_count] + (u[k-1][h_count] - 2 * u[k][h_count]) / (tau * tau));
			}

			u[k+1] = matrix.Solve(v);
		}
		return u;
	}
}
