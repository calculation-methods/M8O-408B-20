#pragma once

#include <functional>
#include <vector>
#include <tuple>
#include <cmath>

#include "common.hpp"

namespace EllipticPDE {
	template <class T>
	using grid_t = std::vector<std::vector<T>>;

	template <class T>
	struct PDE {
		using f_x = std::function<T(T)>;
		using f_y = f_x;
		using f_x_y = std::function<T(T, T)>;

		T a, bx, by, c;
		f_x_y f;
		T x0, x1, y0, y1;
		T alpha_x0, beta_x0;
		f_y gamma_x0;
		T alpha_x1, beta_x1;
		f_y gamma_x1;
		T alpha_y0, beta_y0;
		f_x gamma_y0;
		T alpha_y1, beta_y1;
		f_x gamma_y1;
		f_x_y solution;

		PDE() = default;

		PDE(T a, T bx, T by, T c, f_x_y f, T x0, T x1, T y0, T y1, 
		    T alpha_x0, T beta_x0, f_y gamma_x0, T alpha_x1, T beta_x1, f_y gamma_x1, 
		    T alpha_y0, T beta_y0, f_x gamma_y0, T alpha_y1, T beta_y1, f_x gamma_y1, f_x_y solution) : 
		    a(a), bx(bx), by(by), c(c), f(f), x0(x0), x1(x1), y0(y0), y1(y1), 
		    alpha_x0(alpha_x0), beta_x0(beta_x0), gamma_x0(gamma_x0), alpha_x1(alpha_x1), beta_x1(beta_x1), gamma_x1(gamma_x1), 
		    alpha_y0(alpha_y0), beta_y0(beta_y0), gamma_y0(gamma_y0), alpha_y1(alpha_y1), beta_y1(beta_y1), gamma_y1(gamma_y1), solution(solution) {}

		PDE(T a, T bx, T by, T c, f_x_y f, T x0, T x1, T y0, T y1, 
		    T alpha_x0, T beta_x0, f_y gamma_x0, T alpha_x1, T beta_x1, f_y gamma_x1, 
		    T alpha_y0, T beta_y0, f_x gamma_y0, T alpha_y1, T beta_y1, f_x gamma_y1) : 
		    a(a), bx(bx), by(by), c(c), f(f), x0(x0), x1(x1), y0(y0), y1(y1), 
		    alpha_x0(alpha_x0), beta_x0(beta_x0), gamma_x0(gamma_x0), alpha_x1(alpha_x1), beta_x1(beta_x1), gamma_x1(gamma_x1), 
		    alpha_y0(alpha_y0), beta_y0(beta_y0), gamma_y0(gamma_y0), alpha_y1(alpha_y1), beta_y1(beta_y1), gamma_y1(gamma_y1) {}

		void SetEquation(T a_, T bx_, T by_, T c_, f_x_y f_) {
			a = a_;
			bx = bx_;
			by = by_;
			c = c_;
			f = f_;
		}

		void SetBoundaries(T x0_, T x1_, T y0_, T y1_, T alpha_x0_, T beta_x0_, f_y gamma_x0_, 
		                   T alpha_x1_, T beta_x1_, f_y gamma_x1_, T alpha_y0_, T beta_y0_, f_x gamma_y0_, 
		                   T alpha_y1_, T beta_y1_, f_x gamma_y1_) {
			x0 = x0_; x1 = x1_; y0 = y0_; y1 = y1_;
			alpha_x0 = alpha_x0_; beta_x0 = beta_x0_; gamma_x0 = gamma_x0_;
			alpha_x1 = alpha_x1_; beta_x1 = beta_x1_; gamma_x1 = gamma_x1_;
			alpha_y0 = alpha_y0_; beta_y0 = beta_y0_; gamma_y0 = gamma_y0_;
			alpha_y1 = alpha_y1_; beta_y1 = beta_y1_; gamma_y1 = gamma_y1_;
		}

		void SetSolution(f_x_y solution_) {
			solution = solution_;
		}
	};

	template <class T>
	inline T Relax(T old_value, T new_value, T relax) {
		return old_value + relax * (new_value - old_value);
	}

	template <class T>
	std::tuple<std::vector<T>, std::vector<T>, grid_t<T>, int>
	IterationSolver(const PDE<T>& pde, int nx, int ny, T eps, double relax = 1) {
		std::vector<T> x(nx + 1), y(ny + 1);
		grid_t<T> u(nx + 1, std::vector<T>(ny + 1, 0)), next_u(nx + 1, std::vector<T>(ny + 1, 0));

		T hx = (pde.x1 - pde.x0) / nx;
		T hy = (pde.y1 - pde.y0) / ny;

		for (int i = 0; i <= ny; ++i) {
			y[i] = pde.y0 + i * hy;
		}
		for (int i = 0; i <= nx; ++i) {
			x[i] = pde.x0 + i * hx;
		}

		for (int i = 0; i <= nx; ++i) {
			u[i][0] = pde.gamma_y0(x[i]) / (pde.beta_y0 - pde.alpha_y0 / hy);
			u[i][ny] = pde.gamma_y1(x[i]) / (pde.beta_y1 + pde.alpha_y1 / hy);
		}
		for (int i = 0; i <= ny; ++i) {
			u[0][i] = pde.gamma_x0(y[i]) / (pde.beta_x0 - pde.alpha_x0 / hx);
			u[nx][i] = pde.gamma_x1(y[i]) / (pde.beta_x1 + pde.alpha_x1 / hx);
		}

		T eps_k;
		int iter_count = 0;
		T coeff = (2 * pde.a * (1.0 / hx / hx + 1.0 / hy / hy) + pde.c);
		do {
			eps_k = 0;
			for (int i = 1; i < nx; ++i) {
				for (int j = 1; j < ny; ++j) {
					next_u[i][j] = Relax(u[i][j], 
						                 (pde.a * ((u[i+1][j] + u[i-1][j]) / hx / hx + (u[i][j+1] + u[i][j-1]) / hy / hy) -
					                     (pde.bx * (u[i+1][j] - u[i-1][j]) / hx + pde.by * (u[i][j+1] - u[i][j-1]) / hy) / 2 - pde.f(x[i], y[j])) / coeff, 
					                     relax);

					eps_k = std::max(eps_k, std::abs(next_u[i][j] - u[i][j]));
				}
			}

			for (int i = 1; i < nx; ++i) {
				next_u[i][0] = Relax(u[i][0], (pde.gamma_y0(x[i]) - pde.alpha_y0 / hy * next_u[i][1]) / (pde.beta_y0 - pde.alpha_y0 / hy), relax);
				next_u[i][ny] = Relax(u[i][ny], (pde.gamma_y1(x[i]) + pde.alpha_y1 / hy * next_u[i][ny-1]) / (pde.beta_y1 + pde.alpha_y1 / hy), relax);

				eps_k = std::max(eps_k, std::abs(next_u[i][0] - u[i][0]));
				eps_k = std::max(eps_k, std::abs(next_u[i][ny] - u[i][ny]));
			}
			for (int i = 1; i < ny; ++i) {
				next_u[0][i] = Relax(u[0][i], (pde.gamma_x0(y[i]) - pde.alpha_x0 / hx * next_u[1][i]) / (pde.beta_x0 - pde.alpha_x0 / hx), relax);
				next_u[nx][i] = Relax(u[nx][i], (pde.gamma_x1(y[i]) + pde.alpha_x1 / hx * next_u[nx-1][i]) / (pde.beta_x1 + pde.alpha_x1 / hx), relax);

				eps_k = std::max(eps_k, std::abs(next_u[0][i] - u[0][i]));
				eps_k = std::max(eps_k, std::abs(next_u[nx][i] - u[nx][i]));
			}

			std::swap(next_u, u);

			++iter_count;

		} while (eps_k > eps);

		u[0][0] = (pde.gamma_y0(x[0]) - pde.alpha_y0 / hx * u[0][1]) / (pde.beta_y0 - pde.alpha_y0 / hx);
		u[nx][0] = (pde.gamma_y0(x[nx]) - pde.alpha_y0 / hx * u[nx][1]) / (pde.beta_y0 - pde.alpha_y0 / hx);
		u[0][ny] = (pde.gamma_y1(x[0]) + pde.alpha_y1 / hx * u[0][ny-1]) / (pde.beta_y1 + pde.alpha_y1 / hx);
		u[nx][ny] = (pde.gamma_y1(x[nx]) + pde.alpha_y1 / hx * u[nx][ny-1]) / (pde.beta_y1+  pde.alpha_y1 / hx);

		return {x, y, u, iter_count};
	}

	template <class T>
	std::tuple<std::vector<T>, std::vector<T>, grid_t<T>, int>
	SeidelSolver(const PDE<T>& pde, int nx, int ny, T eps, double relax = 1) {
		std::vector<T> x(nx + 1), y(ny + 1);
		grid_t<T> u(nx + 1, std::vector<T>(ny + 1, 0));

		T hx = (pde.x1 - pde.x0) / nx;
		T hy = (pde.y1 - pde.y0) / ny;

		for (int i = 0; i <= ny; ++i) {
			y[i] = pde.y0 + i * hy;
		}
		for (int i = 0; i <= nx; ++i) {
			x[i] = pde.x0 + i * hx;
		}

		for (int i = 0; i <= nx; ++i) {
			u[i][0] = pde.gamma_y0(x[i]) / (pde.beta_y0 - pde.alpha_y0 / hy);
			u[i][ny] = pde.gamma_y1(x[i]) / (pde.beta_y1 + pde.alpha_y1 / hy);
		}
		for (int i = 0; i <= ny; ++i) {
			u[0][i] = pde.gamma_x0(y[i]) / (pde.beta_x0 - pde.alpha_x0 / hx);
			u[nx][i] = pde.gamma_x1(y[i]) / (pde.beta_x1 + pde.alpha_x1 / hx);
		}

		T eps_k, next_u;
		int iter_count = 0;
		T coeff = (2 * pde.a * (1.0 / hx / hx + 1.0 / hy / hy) + pde.c);
		do {
			eps_k = 0;
			for (int i = 1; i < nx; ++i) {
				for (int j = 1; j < ny; ++j) {
					next_u = Relax(u[i][j], 
						           (pde.a * ((u[i+1][j] + u[i-1][j]) / hx / hx + (u[i][j+1] + u[i][j-1]) / hy / hy) -
					               (pde.bx * (u[i+1][j] - u[i-1][j]) / hx + pde.by * (u[i][j+1] - u[i][j-1]) / hy) / 2 - pde.f(x[i], y[j])) / coeff, 
					               relax);

					eps_k = std::max(eps_k, std::abs(next_u - u[i][j]));
					u[i][j] = next_u;
				}
			}

			for (int i = 1; i < nx; ++i) {
				next_u = Relax(u[i][0], (pde.gamma_y0(x[i]) - pde.alpha_y0 / hy * u[i][1]) / (pde.beta_y0 - pde.alpha_y0 / hy), relax);
				eps_k = std::max(eps_k, std::abs(next_u - u[i][0]));
				u[i][0] = next_u;

				next_u = Relax(u[i][ny], (pde.gamma_y1(x[i]) + pde.alpha_y1 / hy * u[i][ny-1]) / (pde.beta_y1 + pde.alpha_y1 / hy), relax);
				eps_k = std::max(eps_k, std::abs(next_u - u[i][ny]));
				u[i][ny] = next_u;
			}
			for (int i = 1; i < ny; ++i) {
				next_u = Relax(u[0][i], (pde.gamma_x0(y[i]) - pde.alpha_x0 / hx * u[1][i]) / (pde.beta_x0 - pde.alpha_x0 / hx), relax);
				eps_k = std::max(eps_k, std::abs(next_u - u[0][i]));
				u[0][i] = next_u;

				next_u = Relax(u[nx][i], (pde.gamma_x1(y[i]) + pde.alpha_x1 / hx * u[nx-1][i]) / (pde.beta_x1 + pde.alpha_x1 / hx), relax);
				eps_k = std::max(eps_k, std::abs(next_u - u[nx][i]));
				u[nx][i] = next_u;
			}

			++iter_count;

		} while (eps_k > eps);

		u[0][0] = (pde.gamma_y0(x[0]) - pde.alpha_y0 / hx * u[0][1]) / (pde.beta_y0 - pde.alpha_y0 / hx);
		u[nx][0] = (pde.gamma_y0(x[nx]) - pde.alpha_y0 / hx * u[nx][1]) / (pde.beta_y0 - pde.alpha_y0 / hx);
		u[0][ny] = (pde.gamma_y1(x[0]) + pde.alpha_y1 / hx * u[0][ny-1]) / (pde.beta_y1 + pde.alpha_y1 / hx);
		u[nx][ny] = (pde.gamma_y1(x[nx]) + pde.alpha_y1 / hx * u[nx][ny-1]) / (pde.beta_y1+  pde.alpha_y1 / hx);

		return {x, y, u, iter_count};
	}
}