#pragma once

#include <vector>
#include <iostream>
#include <initializer_list>
#include <cmath>
#include <string>

template <class T>
class Polynomial {
private:
	std::vector<T> _data;
	T x0;

	constexpr static double EPS = 1e-9;

public:
	Polynomial(T x0 = 0) : _data(1, T()), x0(x0) { }

	Polynomial(size_t degree, T x0 = 0) : _data(degree+1), x0(x0) { }

	Polynomial(const std::vector<T>& vec, T x0 = 0) : _data(vec), x0(x0) {
		if (vec.size() == 0) {
			_data.push_back(T());
		}
	}

	Polynomial(std::initializer_list<T> list) : _data(list), x0(0) {
		if (list.size() == 0) {
			_data.push_back(T());
		}
	}

	size_t Degree() const {
		return _data.size() - 1;
	}

	size_t Size() const {
		return _data.size();
	}

	T& operator[](size_t index) {
		return _data[index];
	}

	const T& operator[](size_t index) const {
		return _data[index];
	}

	T operator()(T x) {
		T res = 0;
		for (size_t i = Degree(); i > 0; --i) {
			res = (res +  _data[i]) * x;
		}
		res += _data[0];
		return res;
	}

	void Assign(size_t degree, T value = T()) {
		_data.assign(degree + 1, value);
	}

	T Integrate(T a, T b) {
		T res = 0;
		for (size_t i = 0; i < Size(); ++i) {
			res += _data[i] / (i+1) * (std::pow(b, i+1) - std::pow(a, i+1));
		}
		return res;
	}

	template <class U>
	Polynomial<T>& operator+=(const Polynomial<U>& other) {
		size_t len = std::min(_data.size(), other.Size());
		for (size_t i = 0; i < len; ++i) {
			_data[i] += other[i];
		}
		if (_data.size() < other.Size()) {
			for (size_t i = _data.size(); i < other.Size(); ++i) {
				_data.push_back(other[i]);
			}
		}
		return *this;
	}

	template <class U>
	Polynomial<T>& operator*=(const Polynomial<U>& other) {
		std::vector<T> self = _data;
		Assign(Degree() + other.Degree());
		for (size_t i = 0; i < other.Size(); ++i) {
			for (size_t j = 0; j < self.size(); ++j) {
				_data[j + i] += self[j] * other[i];
			}
		}
		return *this;
	}

	template <class U>
	Polynomial<T>& operator*=(const U& value) {
		for (size_t i = 0; i < _data.size(); ++i) {
			_data[i] *= value;
		}
		return *this;
	}

	template <class U>
	friend std::ostream& operator<<(std::ostream&, const Polynomial<U>&);
};

template <class T>
std::ostream& operator<<(std::ostream& os, const Polynomial<T>& p) {
	for (size_t i = p.Degree(); i > 0; --i) {
		if (p[i] < 0) {
			os << "- ";
		} else if (i != p.Degree()) {
			os << "+ ";
		}
		os << std::abs(p[i]) << " * ";
		if (std::abs(p.x0) < Polynomial<T>::EPS) {
			os << "x";
		} else {
			os << "(x - " << p.x0 << ")";
		}
		os << "**" << i << ' ';
	}

	if (p[0] < 0) {
		os << "- ";
	} else if (p.Degree() != 0) {
		os << "+ ";
	}
	os << std::abs(p[0]);
	return os;
}

template<class T, class U>
Polynomial<T> operator*(const Polynomial<T>& a, const Polynomial<U>& b) {
	Polynomial<T> res(a.Degree() + b.Degree());
	for (size_t i = 0; i < a.Size(); ++i) {
		for (size_t j = 0; j < b.Size(); ++j) {
			res[j + i] += b[j] * a[i];
		}
	}
	return res;
}

template<class T, class U>
Polynomial<T> operator+(const Polynomial<T>& a, const Polynomial<U>& b) {
	Polynomial<T> res(std::max(a.Degree(), b.Degree()));
	for (size_t i = 0; i < a.Size(); ++i) {
		res[i] += a[i];
	}
	for (size_t i = 0; i < b.Size(); ++i) {
		res[i] += b[i];
	}
	return res;
}

template<class T, class U>
Polynomial<T> operator+(const Polynomial<T>& a, const U& b) {
	Polynomial<T> res(a);
	res[0] += b;
	return res;
}

template<class T, class U>
Polynomial<T> operator*(const T& a, const Polynomial<U>& b) {
	Polynomial<T> res(b.Degree());
	for (size_t i = 0; i < b.Size(); ++i) {
		res[i] = a * b[i];
	}
	return res;
}