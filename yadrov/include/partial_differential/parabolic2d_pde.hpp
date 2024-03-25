#pragma once

#include <functional>
#include <vector>
#include <tuple>

#include "../linear/tridiagonal_matrix.hpp"
#include "../linear/vector.hpp"
#include "common.hpp"

template <int D, class T>
struct Tensor : public std::vector<Tensor<D - 1, T>> {
	static_assert(D >= 1, "Tensor dimension must be greater than zero");

	template <class... Args>
	Tensor(int size = 0, Args... args) : std::vector<Tensor<D - 1, T>>(size, Tensor<D - 1, T>(args...)) {}
};

template <class T>
struct Tensor<1, T> : public std::vector<T> {
	Tensor(int size = 0, const T& value = T()): std::vector<T>(size, value) {}
};

namespace Parabolic2dPDE {
	template <class T>
	struct PDE {
		using f_x_y = std::function<T(T, T)>;
		using f_x_t = f_x_y;
		using f_y_t = f_x_y;
		using f_x_y_t = std::function<T(T, T, T)>;

		T a, bx, by, c;
		f_x_y_t f;
		T x0, x1, y0, y1;
		f_x_y psi;
		T alpha_x0, beta_x0;
		f_y_t gamma_x0;
		T alpha_x1, beta_x1;
		f_y_t gamma_x1;
		T alpha_y0, beta_y0;
		f_x_t gamma_y0;
		T alpha_y1, beta_y1;
		f_x_t gamma_y1;
		f_x_y_t solution;

		PDE() = default;

		PDE(T a, T bx, T by, T c, f_x_y_t f, T x0, T x1, T y0, T y1, f_x_y psi) : 
		    a(a), bx(bx), by(by), c(c), f(f), x0(x0), x1(x1), y0(y0), y1(y1), psi(psi) {}

		PDE(T a, T bx, T by, T c, f_x_y_t f, T x0, T x1, T y0, T y1, f_x_y psi,
		    T alpha_x0, T beta_x0, f_y_t gamma_x0, T alpha_x1, T beta_x1, f_y_t gamma_x1, 
		    T alpha_y0, T beta_y0, f_x_t gamma_y0, T alpha_y1, T beta_y1, f_x_t gamma_y1, f_x_y_t solution) : 
		    a(a), bx(bx), by(by), c(c), f(f), x0(x0), x1(x1), y0(y0), y1(y1), psi(psi),
		    alpha_x0(alpha_x0), beta_x0(beta_x0), gamma_x0(gamma_x0), alpha_x1(alpha_x1), beta_x1(beta_x1), gamma_x1(gamma_x1), 
		    alpha_y0(alpha_y0), beta_y0(beta_y0), gamma_y0(gamma_y0), alpha_y1(alpha_y1), beta_y1(beta_y1), gamma_y1(gamma_y1), solution(solution) {}

		PDE(T a, T bx, T by, T c, f_x_y_t f, T x0, T x1, T y0, T y1, f_x_y psi,
		    T alpha_x0, T beta_x0, f_y_t gamma_x0, T alpha_x1, T beta_x1, f_y_t gamma_x1, 
		    T alpha_y0, T beta_y0, f_x_t gamma_y0, T alpha_y1, T beta_y1, f_x_t gamma_y1) : 
		    a(a), bx(bx), by(by), c(c), f(f), x0(x0), x1(x1), y0(y0), y1(y1), psi(psi),
		    alpha_x0(alpha_x0), beta_x0(beta_x0), gamma_x0(gamma_x0), alpha_x1(alpha_x1), beta_x1(beta_x1), gamma_x1(gamma_x1), 
		    alpha_y0(alpha_y0), beta_y0(beta_y0), gamma_y0(gamma_y0), alpha_y1(alpha_y1), beta_y1(beta_y1), gamma_y1(gamma_y1) {}

		void SetEquation(T a_, T bx_, T by_, T c_, f_x_y_t f_, T x0_, T x1_, T y0_, T y1_, f_x_y psi_) {
			a = a_;
			bx = bx_;
			by = by_;
			c = c_;
			f = f_;
			x0 = x0_; x1 = x1_; y0 = y0_; y1 = y1_;
			psi = psi_;
		}

		void SetBoundaries(T alpha_x0_, T beta_x0_, f_y_t gamma_x0_, 
		                   T alpha_x1_, T beta_x1_, f_y_t gamma_x1_,
		                   T alpha_y0_, T beta_y0_, f_x_t gamma_y0_, 
		                   T alpha_y1_, T beta_y1_, f_x_t gamma_y1_) {
			alpha_x0 = alpha_x0_; beta_x0 = beta_x0_; gamma_x0 = gamma_x0_;
			alpha_x1 = alpha_x1_; beta_x1 = beta_x1_; gamma_x1 = gamma_x1_;
			alpha_y0 = alpha_y0_; beta_y0 = beta_y0_; gamma_y0 = gamma_y0_;
			alpha_y1 = alpha_y1_; beta_y1 = beta_y1_; gamma_y1 = gamma_y1_;
		}

		void SetSolution(f_x_y_t solution_) {
			solution = solution_;
		}
	};

	template <class T>
	std::tuple<Tensor<3, T>, std::vector<T>, std::vector<T>, std::vector<T>>
	AlternatingDirectionMethod(const PDE<T>& pde, int nx, int ny, int nt, T t_end) {
		std::vector<T> x(nx + 1), y(ny + 1), t(nt + 1);
		Tensor<3, T> u(nt + 1, nx + 1, ny + 1);
		Tensor<2, T> p(nx + 1, ny + 1);

		T hx = (pde.x1 - pde.x0) / nx;
		T hy = (pde.y1 - pde.y0) / ny;
		T tau = t_end / nt;

		for (int i = 0; i <= nx; ++i) {
			x[i] = pde.x0 + i * hx;
		}
		for (int i = 0; i <= ny; ++i) {
			y[i] = pde.y0 + i * hy;
		}
		for (int i = 0; i <= nt; ++i) {
			t[i] = i * tau;
		}

		for (size_t i = 0; i < x.size(); ++i) {
			for (size_t j = 0; j < y.size(); ++j) {
				u[0][i][j] = pde.psi(x[i], y[j]);
			}
		}

		Vector<T> vx(nx+1), vy(ny+1);
		TDMatrix<T> mx(nx+1), my(ny+1);

		T alpha = (- pde.a / hx + pde.bx / 2) / hx,
		  beta = 2 * pde.a / hx / hx + 2 / tau - pde.c,
		  gamma = (- pde.a / hx - pde.bx / 2) / hx;
		for (int i = 1; i < nx; ++i) {
			mx.a[i] = alpha;
			mx.b[i] = beta;
			mx.c[i] = gamma;
		}

		alpha = (- pde.a / hy + pde.by / 2) / hy,
		beta = 2 * pde.a / hy / hy + 2 / tau - pde.c,
		gamma = (- pde.a / hy - pde.by / 2) / hy;
		for (int j = 1; j < ny; ++j) {
			my.a[j] = alpha;
			my.b[j] = beta;
			my.c[j] = gamma;
		}

		for (int k = 0; k < nt; ++k) {
			for (int j = 1; j < ny; ++j) {
				for (int i = 1; i < nx; ++i) {
					vx[i] = u[k][i][j-1] * (pde.a / hy - pde.by / 2) / hy + 
					        2 * u[k][i][j] * (1 / tau - pde.a / hy / hy) + 
					        u[k][i][j+1] * (pde.a / hy + pde.by / 2) / hy +
					        pde.f(x[i], y[j], t[k] + tau / 2);
				}
				vx[0] = pde.gamma_x0(y[j], t[k] + tau / 2);
				vx[nx] = pde.gamma_x1(y[j], t[k] + tau / 2);

				mx.b[0] = -pde.alpha_x0 / hx + pde.beta_x0;
				mx.c[0] = pde.alpha_x0 / hx;

				mx.a[nx] = -pde.alpha_x1 / hx;
				mx.b[nx] = pde.alpha_x1 / hx + pde.beta_x1;

				vx = mx.Solve(vx);

				for (int i = 0; i <= nx; ++i) {
					p[i][j] = vx[i];
				}
			}

			for (int i = 0; i <= nx; ++i) {
				p[i][0] = (pde.gamma_y0(x[i], t[k] + tau/2) - pde.alpha_y0 / hy * p[i][1]) / (-pde.alpha_y0 / hy + pde.beta_y0);
				p[i][ny] = (pde.gamma_y1(x[i], t[k] + tau/2) + pde.alpha_y1 / hy * p[i][ny-1]) / (pde.alpha_y1 / hy + pde.beta_y1);
			}

			for (int i = 1; i < nx; ++i) {
				for (int j = 1; j < ny; ++j) {
					vy[j] = p[i-1][j] * (pde.a / hx - pde.bx / 2) / hx + 
					        2 * p[i][j] * (1 / tau - pde.a / hx / hx) + 
					        p[i+1][j] * (pde.a / hx + pde.bx / 2) / hx +
					        pde.f(x[i], y[j], t[k+1]);
				}
				vy[0] = pde.gamma_y0(x[i], t[k+1]);
				vy[ny] = pde.gamma_y1(x[i], t[k+1]);

				my.b[0] = -pde.alpha_y0 / hy + pde.beta_y0;
				my.c[0] = pde.alpha_y0 / hy;

				my.a[ny] = -pde.alpha_y1 / hy;
				my.b[ny] = pde.alpha_y1 / hy + pde.beta_y1;

				vy = my.Solve(vy);

				for (int j = 0; j <= ny; ++j) {
					u[k+1][i][j] = vy[j];
				}
			}

			for (int j = 0; j <= ny; ++j) {
				u[k+1][0][j] = (pde.gamma_x0(y[j], t[k+1]) - pde.alpha_x0 / hx * u[k+1][1][j]) / (-pde.alpha_x0 / hx + pde.beta_x0);
				u[k+1][nx][j] = (pde.gamma_x1(y[j], t[k+1]) + pde.alpha_x1 / hx * u[k+1][nx-1][j]) / (pde.alpha_x1 / hx + pde.beta_x1);
			}
		}

		return {u, x, y, t};
	}

	template <class T>
	std::tuple<Tensor<3, T>, std::vector<T>, std::vector<T>, std::vector<T>>
	FractionalStepMethod(const PDE<T>& pde, int nx, int ny, int nt, T t_end) {
		std::vector<T> x(nx + 1), y(ny + 1), t(nt + 1);
		Tensor<3, T> u(nt + 1, nx + 1, ny + 1);
		Tensor<2, T> p(nx + 1, ny + 1);

		T hx = (pde.x1 - pde.x0) / nx;
		T hy = (pde.y1 - pde.y0) / ny;
		T tau = t_end / nt;

		for (int i = 0; i <= nx; ++i) {
			x[i] = pde.x0 + i * hx;
		}
		for (int i = 0; i <= ny; ++i) {
			y[i] = pde.y0 + i * hy;
		}
		for (int i = 0; i <= nt; ++i) {
			t[i] = i * tau;
		}

		for (size_t i = 0; i < x.size(); ++i) {
			for (size_t j = 0; j < y.size(); ++j) {
				u[0][i][j] = pde.psi(x[i], y[j]);
			}
		}

		Vector<T> vx(nx+1), vy(ny+1);
		TDMatrix<T> mx(nx+1), my(ny+1);

		T alpha = (- pde.a / hx + pde.bx / 2) / hx,
		  beta = 2 * pde.a / hx / hx + 1 / tau - pde.c,
		  gamma = (- pde.a / hx - pde.bx / 2) / hx;
		for (int i = 1; i < nx; ++i) {
			mx.a[i] = alpha;
			mx.b[i] = beta;
			mx.c[i] = gamma;
		}

		alpha = (- pde.a / hy + pde.by / 2) / hy,
		beta = 2 * pde.a / hy / hy + 1 / tau - pde.c,
		gamma = (- pde.a / hy - pde.by / 2) / hy;
		for (int j = 1; j < ny; ++j) {
			my.a[j] = alpha;
			my.b[j] = beta;
			my.c[j] = gamma;
		}

		for (int k = 0; k < nt; ++k) {
			for (int j = 1; j < ny; ++j) {
				for (int i = 1; i < nx; ++i) {
					vx[i] = u[k][i][j] / tau + pde.f(x[i], y[j], t[k]) / 2;
				}
				vx[0] = pde.gamma_x0(y[j], t[k+1]);
				vx[nx] = pde.gamma_x1(y[j], t[k+1]);

				mx.b[0] = -pde.alpha_x0 / hx + pde.beta_x0;
				mx.c[0] = pde.alpha_x0 / hx;

				mx.a[nx] = -pde.alpha_x1 / hx;
				mx.b[nx] = pde.alpha_x1 / hx + pde.beta_x1;

				vx = mx.Solve(vx);

				for (int i = 0; i <= nx; ++i) {
					p[i][j] = vx[i];
				}
			}

			for (int i = 0; i <= nx; ++i) {
				p[i][0] = (pde.gamma_y0(x[i], t[k+1]) - pde.alpha_y0 / hy * p[i][1]) / (-pde.alpha_y0 / hy + pde.beta_y0);
				p[i][ny] = (pde.gamma_y1(x[i], t[k+1]) + pde.alpha_y1 / hy * p[i][ny-1]) / (pde.alpha_y1 / hy + pde.beta_y1);
			}

			for (int i = 1; i < nx; ++i) {
				for (int j = 1; j < ny; ++j) {
					vy[j] = p[i][j] / tau + pde.f(x[i], y[j], t[k+1]) / 2;
				}
				vy[0] = pde.gamma_y0(x[i], t[k+1]);
				vy[ny] = pde.gamma_y1(x[i], t[k+1]);

				my.b[0] = -pde.alpha_y0 / hy + pde.beta_y0;
				my.c[0] = pde.alpha_y0 / hy;

				my.a[ny] = -pde.alpha_y1 / hy;
				my.b[ny] = pde.alpha_y1 / hy + pde.beta_y1;

				vy = my.Solve(vy);

				for (int j = 0; j <= ny; ++j) {
					u[k+1][i][j] = vy[j];
				}
			}

			for (int j = 0; j <= ny; ++j) {
				u[k+1][0][j] = (pde.gamma_x0(y[j], t[k+1]) - pde.alpha_x0 / hx * u[k+1][1][j]) / (-pde.alpha_x0 / hx + pde.beta_x0);
				u[k+1][nx][j] = (pde.gamma_x1(y[j], t[k+1]) + pde.alpha_x1 / hx * u[k+1][nx-1][j]) / (pde.alpha_x1 / hx + pde.beta_x1);
			}
		}

		return {u, x, y, t};
	}
}