#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

const double pi = 3.14159265359;

double a, b, c, d, sigma, N, timeSlice, K, h, tau, mu;

double phi0(double t) {
	return exp(-t);
}

double phi1(double t) {
	return -phi0(t);
}

double psi0(double x) {
	return sin(x);
}

double psi1(double x) {
	return -sin(x);
}

double func(double x, double t) {
	return -cos(x) * exp(-t);
}

double analSol(double x, double t) {
	return exp(-t) * sin(x);
}

double deriv(double (*f)(double), double x) {
	double dx = 0.000001;
	return (f(x + dx) - f(x)) / dx;
}
double deriv2(double (*innerFunc)(double x1), double x) {
	double dx = 0.000001;
	return (deriv(innerFunc, x + dx) - deriv(innerFunc, x)) / dx;
}

void init1Order(vector<vector<double>>& U) {
	for (int j = 0; j < N; ++j) {
		U[1][j] = psi0(j * h) + psi1(j * h) * tau;
	}
}

void implicitInit1Order(vector<vector<double>>& U) {
	for (int j = 0; j < N; ++j) {
		U[1][j] = psi0(j * h) + psi1(j * h) * tau;
	}
}

void init2Order(vector<vector<double>>& U) {
	for (int j = 0; j < N; ++j) {
		U[1][j] = psi0(j * h) + psi1(j * h) * (tau - d * (tau * tau) / 2) + (a *  deriv2(&psi0, j * h) + b * deriv(&psi0, j * h) +
			c * psi0(j * h) + func(j * h, tau)) * (tau * tau) / 2;
	} 
}

void explicitBound2Points1Order(vector<vector<double>>& U) {
	for (int k = 1; k < K - 1; ++k) {
		U[k + 1][0] = (-1.0 / h * U[k + 1][1] + phi0(tau * (k + 1))) / (0.0 - 1.0 / h);
		U[k + 1][N - 1] = (1.0 / h * U[k + 1][N - 2] + phi1(tau * (k + 1))) / (0.0 + 1.0 / h);
	}
}
void explicitBound3Points2Order(vector<vector<double>>& U) {
	for (int k = 1; k < K - 1; ++k) {
		U[k + 1][0] = (-1.0 / h / 2 * (4 * U[k + 1][1] - U[k + 1][2]) + phi0(tau * (k + 1))) /
			(0.0 - 3 * 1.0 / h / 2);
		U[k + 1][N - 1] = (-1.0 / h / 2 * (U[k + 1][N - 3] - 4 * U[k + 1][N - 2]) + phi1(tau * (k + 1))) /
			(0.0 + 3 * 1.0 / h / 2);
	}
}

void explicitBound2Points2Order(vector<vector<double>>& U) {
	double b0 = 2.0 * a / h + h / (tau * tau) - h * c - 0.0 / 1.0 * (2 * a - b * h) + d * h / 2 / tau;
	double c0 = -2.0 * a / h;
	double bN = 2.0 * a / h + h / (tau * tau) - h * c + 0.0 / 1.0 * (2 * a + b * h) + d * h / 2 / tau;
	double aN = -2.0 * a / h;
	for (int k = 1; k < K - 1; ++k) {
		double d0 = 2 * h / (tau * tau) * U[k][0] + (-h / (tau * tau) + d * h / 2 / tau) * U[k - 1][0] -
			phi0(tau * (k + 1)) * (2 * a - b * h) / 1.0 + h * func(0, tau * (k + 1));
		double dN = 2 * h / (tau * tau) * U[k][N - 1] + (-h / (tau * tau) + d * h / 2 / tau) * U[k - 1][N - 1] +
			phi1(tau * (k + 1)) * (2 * a + b * h) / 1.0 - h * func(pi, tau * (k + 1));
		U[k + 1][0] = (d0 - c0 * U[k + 1][1]) / b0;
		U[k + 1][N - 1] = (dN - aN * U[k + 1][N - 2]) / bN;
	}
}

vector<vector<double>> explicit2Points1Bound1Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
    init1Order(U);
	for (int k = 1; k < K - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			U[k + 1][j] = ((sigma + mu) * U[k][j + 1] + (-2 * sigma + 2 + c * tau * tau) * U[k][j] + (sigma - mu) * U[k][j - 1]
				+ (-1 + d * tau / 2) * U[k - 1][j] + tau * tau * func(j * h, k * tau)) / (1 + d * tau / 2);
		}
	}
    explicitBound2Points1Order(U);
	return U;
}

vector<vector<double>> explicit2Points1Bound2Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
    init2Order(U);
	for (int k = 1; k < K - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			U[k + 1][j] = ((sigma + mu) * U[k][j + 1] + (-2 * sigma + 2 + c * tau * tau) * U[k][j] + (sigma - mu) * U[k][j - 1]
				+ (-1 + d * tau / 2) * U[k - 1][j] + tau * tau * func(j * h, k * tau)) / (1 + d * tau / 2);
		}
	}
	explicitBound3Points2Order(U);
	return U;
}

vector<vector<double>> explicit3Points2Bound1Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
    init1Order(U);
	for (int k = 1; k < K - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			U[k + 1][j] = ((sigma + mu) * U[k][j + 1] + (-2 * sigma + 2 + c * tau * tau) * U[k][j] + (sigma - mu) * U[k][j - 1]
				+ (-1 + d * tau / 2) * U[k - 1][j] + tau * tau * func(j * h, k * tau)) / (1 + d * tau / 2);
		}
	}
    explicitBound3Points2Order(U);
	return U;
}

vector<vector<double>> explicit3Points2Bound2Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
    init2Order(U);
	for (int k = 1; k < K - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			U[k + 1][j] = ((sigma + mu) * U[k][j + 1] + (-2 * sigma + 2 + c * tau * tau) * U[k][j] + (sigma - mu) * U[k][j - 1]
				+ (-1 + d * tau / 2) * U[k - 1][j] + tau * tau * func(j * h, k * tau)) / (1 + d * tau / 2);
		}
	}
    explicitBound3Points2Order(U);
	return U;
}

vector<vector<double>> explicit2Points2Bound1Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
    init1Order(U);
	for (int k = 1; k < K - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			U[k + 1][j] = ((sigma + mu) * U[k][j + 1] + (-2 * sigma + 2 + c * tau * tau) * U[k][j] + (sigma - mu) * U[k][j - 1]
				+ (-1 + d * tau / 2) * U[k - 1][j] + tau * tau * func(j * h, k * tau)) / (1 + d * tau / 2);
		}
	}
    explicitBound2Points2Order(U);
	return U;
}

vector<vector<double>> explicit2Points2Bound2Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
    init2Order(U);
	for (int k = 1; k < K - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			U[k + 1][j] = ((sigma + mu) * U[k][j + 1] + (-2 * sigma + 2 + c * tau * tau) * U[k][j] + (sigma - mu) * U[k][j - 1]
				+ (-1 + d * tau / 2) * U[k - 1][j] + tau * tau * func(j * h, k * tau)) / (1 + d * tau / 2);
		}
	}
    explicitBound2Points2Order(U);
	return U;
}

double tridiagonalAlgo(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d,
	vector<double>& x, int step, double prevP, double prevQ) {
	if (step == a.size()) {
		return 0;
	}
	double p = ((-1.0) * c[step]) / (b[step] + a[step] * prevP);
	double q = (d[step] - a[step] * prevQ) / (b[step] + a[step] * prevP);
	x[step] = p * tridiagonalAlgo(a, b, c, d, x, step + 1, p, q) + q;
	return x[step];
}

void implicitBound2Points1Order(vector<vector<double>>& U,
	vector<double>& lower, vector<double>& main, vector<double>& upper,
	vector<double>& coeffs, bool getA, int k) {
	double aCoeff = a / (h * h);
	double bCoeff = b / 2.0 / h;
	double tCoeff = 1.0 / (tau * tau);
	double dCoeff = d / 2.0 / tau;
	double aj = bCoeff - aCoeff;
	double bj = tCoeff + dCoeff + 2.0 * aCoeff;
	double cj = -bCoeff - aCoeff;
	if (getA) {
		double b0 = 0.0 - 1.0 / h;
		double c0 = 1.0 / h;
		double aN = -1.0 / h;
		double bN = 0.0 + 1.0 / h;
		lower[0] = 0;
		main[0] = b0;
		upper[0] = c0;
		for (int j = 1; j < N - 1; ++j) {
			lower[j] = aj;
			main[j] = bj;
			upper[j] = cj;
		}
		lower[N - 1] = aN;
		main[N - 1] = bN;
		upper[N - 1] = 0.0;
	}
	else {
		coeffs[0] = phi0(tau * (k + 1));
		coeffs[N - 1] = phi1(tau * (k + 1));
		for (int j = 1; j < N - 1; ++j) {
			coeffs[j] = U[k - 1][j] * (-tCoeff + dCoeff) + U[k][j] * (2.0 * tCoeff + c) + func(j * h, k * tau);
		}
	}
}

void implicitBound3Points2Order(vector<vector<double>>& U,
	vector<double>& lower, vector<double>& main, vector<double>& upper,
	vector<double>& coeffs, bool getA, int k) {
	double aCoeff = a / (h * h);
	double bCoeff = b / 2.0 / h;
	double tCoeff = 1.0 / (tau * tau);
	double dCoeff = d / 2.0 / tau;

	double aj = bCoeff - aCoeff;
	double bj = tCoeff + dCoeff + 2.0 * aCoeff;
	double cj = -bCoeff - aCoeff;
	if (getA) {
		double a0 = 0.0 - 3.0 * 1.0 / h / 2.0;
		double b0 = 2.0 * 1.0 / h;
		double c0 = -1.0 / h / 2.0;
		double aN = 1.0 / h / 2.0;
		double bN = -2.0 * 1.0 / h;
		double cN = 0.0 + 3.0 * 1.0 / h / 2.0;

		lower[0] = a0;
		main[0] = b0;
		upper[0] = c0;
		for (int j = 1; j < N - 1; ++j) {
			lower[j] = aj;
			main[j] = bj;
			upper[j] = cj;
		}
		lower[N - 1] = aN;
		main[N - 1] = bN;
		upper[N - 1] = cN;
	}
	else {
		coeffs[0] = phi0(tau * (k + 1));
		coeffs[N - 1] = phi1(tau * (k + 1));
		for (int j = 1; j < N - 1; ++j) {
			coeffs[j] = U[k - 1][j] * (-tCoeff + dCoeff) + U[k][j] * (2 * tCoeff + c) + func(j * h, k * tau);
		}
	}
}

void implicitBound2Points2Order(vector<vector<double>>& U,
	vector<double>& lower, vector<double>& main, vector<double>& upper,
	vector<double>& coeffs, bool getA, int k) {
	double aCoeff = a / (h * h);
	double bCoeff = b / 2.0 / h;
	double tCoeff = 1.0 / (tau * tau);
	double dCoeff = d / 2.0 / tau;

	double aj = bCoeff - aCoeff;
	double bj = tCoeff + dCoeff + 2.0 * aCoeff;
	double cj = -bCoeff - aCoeff;
	if (getA) { 
		double b0 = 2.0 * a / h + h / (tau * tau) - h * c - 0.0 / 1.0 * (2.0 * a - b * h) + d * h / 2 / tau;
		double c0 = -2.0 * a / h;
		double aN = -2.0 * a / h;
		double bN = 2.0 * a / h + h / (tau * tau) - h * c + 0.0 / 1.0 * (2.0 * a + b * h) + d * h / 2 / tau;

		lower[0] = 0;
		main[0] = b0;
		upper[0] = c0;
		for (int j = 1; j < N - 1; ++j) {
			lower[j] = aj;
			main[j] = bj;
			upper[j] = cj;
		}
		lower[N - 1] = aN;
		main[N - 1] = bN;
		upper[N - 1] = 0;
	} else {
		coeffs[0] = 2 * h / (tau * tau) * U[k][0] + (-h / (tau * tau) + d * h / 2 / tau) * U[k - 1][0] - phi0(tau * (k + 1)) * (
			2 * a - b * h) / 1.0 + h * func(0, tau * (k + 1));
		coeffs[N - 1] = 2 * h / (tau * tau) * U[k][N - 1] + (-h / (tau * tau) + d * h / 2 / tau) * U[k - 1][N - 1] + phi1(
			tau * (k + 1)) * (2 * a + b * h) / 1.0 - h * func(pi, tau * (k + 1));
		for (int j = 1; j < N - 1; ++j) {
			coeffs[j] = U[k - 1][j] * (-tCoeff + dCoeff) + U[k][j] * (2 * tCoeff + c) + func(j * h, k * tau);
		}
	}
}

vector<vector<double>> implicit2Points1Bound1Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
	init1Order(U);
	vector<double> l(N), m(N), u(N), coeffs(N);
	implicitBound2Points1Order(U, l, m, u, coeffs, true, 0);
	for (int k = 1; k < K - 1; ++k) {
		implicitBound2Points1Order(U, l, m, u, coeffs, false, k);
		tridiagonalAlgo(l, m, u, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}

vector<vector<double>> implicit2Points1Bound2Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
	init2Order(U);
	vector<double> l(N), m(N), u(N), coeffs(N);
	implicitBound2Points1Order(U, l, m, u, coeffs, true, 0);
	for (int k = 1; k < K - 1; ++k) {
		implicitBound2Points1Order(U, l, m, u, coeffs, false, k);
		tridiagonalAlgo(l, m, u, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}

vector<vector<double>> implicit3Points2Bound1Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
	init1Order(U);
	vector<double> l(N), m(N), u(N), coeffs(N);
	implicitBound3Points2Order(U, l, m, u, coeffs, true, 0);
	for (int k = 1; k < K - 1; ++k) {
		implicitBound3Points2Order(U, l, m, u, coeffs, false, k);
		tridiagonalAlgo(l, m, u, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}

vector<vector<double>> implicit3Points2Bound2Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
	init2Order(U);
	vector<double> l(N), m(N), u(N), coeffs(N);
	implicitBound3Points2Order(U, l, m, u, coeffs, true, 0);
	for (int k = 1; k < K - 1; ++k) {
		implicitBound3Points2Order(U, l, m, u, coeffs, false, k);
		tridiagonalAlgo(l, m, u, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}

vector<vector<double>> implicit2Points2Bound1Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
	init1Order(U);
	vector<double> l(N), m(N), u(N), coeffs(N);
	implicitBound2Points2Order(U, l, m, u, coeffs, true, 0);
	for (int k = 1; k < K - 1; ++k) {
		implicitBound2Points2Order(U, l, m, u, coeffs, false, k);
		tridiagonalAlgo(l, m, u, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}

vector<vector<double>> implicit2Points2Bound2Init() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
	init2Order(U);
	vector<double> l(N), m(N), u(N), coeffs(N);
	implicitBound2Points2Order(U, l, m, u, coeffs, true, 0);
	for (int k = 1; k < K - 1; ++k) {
		implicitBound2Points2Order(U, l, m, u, coeffs, false, k);
		tridiagonalAlgo(l, m, u, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}

vector<double> getError(vector<vector<double>>& U) {
	vector<double> error(K, 0.0);
	for (int k = 0; k < K; ++k) {
		for (int j = 0; j < N; ++j) {
			double curError = abs(U[k][j] - analSol(h * j, tau * k));
			if (curError > error[k]) {
				error[k] = curError;
			}
		}
	}
	return error;
}

void writeVectorToFile(vector<vector<double>>& U, const string& fileName) {
	ofstream file(fileName);
	vector<double> error = getError(U);
	if (file.is_open()) {
		file << a << " " << b << " " << c << " " << d << " " << sigma << " " << N << " " << timeSlice << endl;
		for (int k = 0; k < K; ++k) {
			for (int j = 0; j < N; ++j) {
				file << U[k][j] << " ";
			}
			file << endl;
		}
		for (int k = 0; k < K; ++k) {
			file << error[k] << ' ';
		}
		file.close();
		cout << "Data written to file: " << fileName << endl;
	} else {
		cerr << "Unable to open file: " << fileName << endl;
	}
}

int main() {
	double pi = 3.14159265359;
	a = 1.0;
	b = 1.0;
	c = -1.0;
    d = 3.0;
	sigma = 0.5;
	N = 100;
	timeSlice = 20;
	K = 500;
	h = pi / N;
	tau = sqrt(sigma * h * h / a);
	mu = b * tau * tau / 2.0 / h;
	vector<vector<double>> U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U11, U12;
	vector<string> fileNames{ "explicit2Points1Bound1Init.txt", 
                                "explicit2Points1Bound2Init.txt", 
                                "explicit3Points2Bound1Init.txt", 
                                "explicit3Points2Bound2Init.txt", 
                                "explicit2Points2Bound1Init.txt", 
                                "explicit2Points2Bound2Init.txt",
                                "implicit2Points1Bound1Init.txt",
                                "implicit2Points1Bound2Init.txt",
                                "implicit3Points2Bound1Init.txt",
                                "implicit3Points2Bound2Init.txt",
                                "implicit2Points2Bound1Init.txt",
                                "implicit2Points2Bound2Init.txt"
                            };
	U1 = explicit2Points1Bound1Init();
	U2 = explicit2Points1Bound2Init();
	U3 = explicit3Points2Bound1Init();
	U4 = explicit3Points2Bound2Init();
	U5 = explicit2Points2Bound1Init();
	U6 = explicit2Points2Bound2Init();
    U7 = implicit2Points1Bound1Init();
    U8 = implicit2Points1Bound2Init();
    U9 = implicit3Points2Bound1Init();
    U10 = implicit3Points2Bound2Init();
    U11 = implicit2Points2Bound1Init();
    U12 = implicit2Points2Bound2Init();

    writeVectorToFile(U1, fileNames[0]);
	writeVectorToFile(U2, fileNames[1]);
	writeVectorToFile(U3, fileNames[2]);
	writeVectorToFile(U4, fileNames[3]);
	writeVectorToFile(U5, fileNames[4]);
	writeVectorToFile(U6, fileNames[5]);
    writeVectorToFile(U7, fileNames[6]);
    writeVectorToFile(U8, fileNames[7]);
    writeVectorToFile(U9, fileNames[8]);
    writeVectorToFile(U10, fileNames[9]);
    writeVectorToFile(U11, fileNames[10]);
    writeVectorToFile(U12, fileNames[11]);
}
