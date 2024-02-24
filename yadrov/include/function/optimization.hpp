// Оффтоп: методы оптимизации

#pragma once

#include <functional>
#include <random>
#include <cmath>

#include "../linear/vector.hpp"

template <class T, class Compare = std::less<T>>
Vector<T> NelderMeadMethod(const std::function<T(Vector<T>)>& f, std::vector<Vector<T>> simplex, T eps, const std::function<bool(Vector<T>)>& region = [](Vector<T>){ return true; }, T alpha = 1, T beta = 0.5, T gamma = 2) {
	for (const Vector<T>& point: simplex) {
		if (point.Size() != simplex.size() - 1) {
			throw std::runtime_error("Dimension mismatch");
		}
	}

	Compare comp;
	T cur_eps = eps + 1;
	size_t best_i, worst_i;
	T best, worst, center_value, reflect_value, value, l, r, cur;
	Vector<T> center(simplex[0].Size()), reflect(center.Size()), point(center.Size());
	size_t iter_count = 0;

	while (cur_eps > eps) {
		best_i = worst_i = 0;
		best = worst = f(simplex[0]);
		for (size_t i = 0; i < simplex.size(); ++i) {
			value = f(simplex[i]);
			if (comp(value, best)) {
				best = value;
				best_i = i;
			}
			if (comp(worst, value)) {
				worst = value;
				worst_i = i;
			}
		}

		center.Fill();
		for (size_t i = 0; i < simplex.size(); ++i) {
			if (i != worst_i) {
				center += simplex[i];
			}
		}
		center /= (simplex.size() - 1);
		center_value = f(center);

		reflect = center + (center - simplex[worst_i]) * alpha;
		if (!region(reflect)) {
			l = 0;
			r = alpha;
			while (r - l > eps) {
				cur = (r + l) / 2;
				if (region(center + (center - simplex[worst_i]) * cur)) {
					l = cur;
				} else {
					r = cur;
				}
			}
			reflect = center + (center - simplex[worst_i]) * l;
		}
		reflect_value = f(reflect);

		if (comp(reflect_value, best)) {
			point = center + (reflect - center) * gamma;
			if (!region(point)) {
				l = 0;
				r = gamma;
				while (r - l > eps) {
					cur = (r + l) / 2;
					if (region(center + (reflect - center) * cur)) {
						l = cur;
					} else {
						r = cur;
					}
				}
				point = center + (reflect - center) * l;
			}

			value = f(point);
			if (comp(value, reflect_value)) {
				simplex[worst_i] = point;
			} else {
				simplex[worst_i] = reflect;
			}

		} else if (comp(reflect_value, worst)) {
			point = center + (simplex[worst_i] - center) * beta;
			value = f(point);
			if (comp(value, reflect_value)) {
				simplex[worst_i] = point;
			} else {
				simplex[worst_i] = reflect;
			}

		} else {
			for (size_t i = 0; i < simplex.size(); ++i) {
				if (i != best_i) {
					simplex[i] = simplex[best_i] + (simplex[i] - simplex[best_i]) * 0.5;
				}
			}
		}

		cur_eps = 0;
		for (size_t i = 0; i < simplex.size(); ++i) {
			cur_eps += std::pow(f(simplex[i]) - center_value, 2);
		}
		cur_eps = std::sqrt(cur_eps / (simplex.size() + 1));
		++iter_count;
	}

	for (size_t i = 0; i < simplex.size(); ++i) {
		value = f(simplex[i]);
		if (comp(value, best)) {
			best = value;
			best_i = i;
		}
	}

	return simplex[best_i];
}

template <class T, class Compare = std::less<T>>
Vector<T> RandomSearch(const std::function<T(Vector<T>)>& f, const Vector<T>& x0, T eps, T r, size_t count, size_t max_test_count, const std::function<bool(Vector<T>)>& region = [](Vector<T>){ return true; }) {
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<T> rand(-1.0, 1.0);
	Compare comp;

	Vector<T> x = x0, next_x(x.Size(), 1), test_x(x.Size()), rand_v(x.Size());
	T min_value, value;
	bool finded;
	size_t iter_count = 0, test_count = 0;
	while (r > eps) {
		min_value = f(x);
		finded = false;

		for (size_t i = 0; i < count; ++i) {
			for (size_t j = 0; j < rand_v.Size(); ++j) {
				rand_v[j] = rand(gen);
			}
			rand_v *= r / rand_v.Norm(2);

			if (!region(x + rand_v)) {
				continue;
			}

			value = f(x + rand_v);
			if (comp(value, min_value)) {
				finded = true;

				T l = 0, r = 10, a, b;
				while (r - l > eps) {
					a = l + (r - l) / 3;
					b = l + 2 * (r - l) / 3;
					if (comp(f(x + rand_v * a), f(x + rand_v * b)) || !region(x + rand_v * b)) {
						r = b;
					} else {
						l = a;
					}
				}

				next_x = rand_v * l;
				min_value = f(x + next_x);
			}
		}

		if (finded) {
			test_count = 0;
			x = x + next_x;
		} else if (test_count < max_test_count) {
			++test_count;
		} else {
			test_count = 0;
			r /= 2;
		}
		++iter_count;
	}
	return x;
}

template <class T, class Compare = std::less<T>>
T TrivialSearch(const Vector<T>& a, const Vector<T>& b, const std::function<T(Vector<T>)>& f, int steps, T min_max, double eps, size_t i = 0, Vector<T> x = {0, 0}) {
	Compare comp;
	T d = (b[i] - a[i]) / steps;
	T max_min = min_max, cur;
	for (T s = a[i]; s < b[i] + eps; s += d) {
		x[i] = s;
		if (i == a.Size() - 1) {
			cur = f(x);
			if (comp(cur, max_min)) {
				max_min = cur;
			}
		} else {
			cur = TrivialSearch<T, Compare>(a, b, f, steps, min_max, eps, i+1, x);
			if (comp(cur, max_min)) {
				max_min = cur;
			}
		}
	}
	return max_min;
}