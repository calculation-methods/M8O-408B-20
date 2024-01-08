#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

double a, b, c, sigma, N, timeSlice, K, h, tau, mu;

double phi0(double t) {
	return exp((c - a) * t) * (cos(b * t) + sin(b * t));
}

double phi1(double t) {
	return -phi0(t);
}

double analSol(double x, double t) {
	return exp((c - a) * t) * sin(x + b * t);
}

vector<vector<double>> explicit2Points1Order() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = sin(j * h);
	}
	for (int k = 0; k < K - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			U[k + 1][j] = (sigma + mu) * U[k][j + 1] + (1 - 2 * sigma + tau * c) * U[k][j] + (sigma - mu) * U[k][j - 1];
		}
	}
	for (int k = 0; k < K - 1; ++k) {
		U[k + 1][0] = (-1.0 / h * U[k + 1][1] + phi0(tau * (k + 1))) / (1.0 - 1.0 / h);
		U[k + 1][N - 1] = (1.0 / h * U[k + 1][N - 2] + phi1(tau * (k + 1))) / (1.0 + 1.0 / h);
	}
	return U;
}

vector<vector<double>> explicit3Points2Order() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = sin(j * h);
	}
	for (int k = 0; k < K - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			U[k + 1][j] = (sigma + mu) * U[k][j + 1] + (1 - 2 * sigma + tau * c) * U[k][j] + (sigma - mu) * U[k][j - 1];
		}
	}
	for (int k = 0; k < K - 1; ++k) {
		U[k + 1][0] = (-1.0 / h / 2 * (4 * U[k + 1][1] - U[k + 1][2]) + phi0(tau * (k + 1))) /
					(1.0 - 3 * 1.0 / h / 2);
		U[k + 1][N - 1] = (-1.0 / h / 2 * (U[k + 1][N - 3] - 4 * U[k + 1][N - 2]) + phi1(tau * (k + 1))) /
						(1.0 + 3 * 1.0 / h / 2);
	}
	return U;
}

vector<vector<double>> explicit2Points2Order() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = sin(j * h);
	}
	for (int k = 0; k < K - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			U[k + 1][j] = (sigma + mu) * U[k][j + 1] + (1 - 2 * sigma + tau * c) * U[k][j] + (sigma - mu) * U[k][j - 1];
		}
	}
	double b0 = 2.0 * a / h + h / tau - h * c - 1.0 / 1.0 * (2 * a - b * h);
	double c0 = -2.0 * a / h;
	double bN = 2.0 * a / h + h / tau - h * c + 1.0 / 1.0 * (2 * a + b * h);
	double aN = -2.0 * a / h;
	for (int k = 0; k < K - 1; ++k) {
		double d0 = h / tau * U[k][0] - phi0(tau * (k + 1)) * (2 * a - b * h) / 1.0;
		double dN = h / tau * U[k][N - 1] + phi1(tau * (k + 1)) * (2 * a + b * h) / 1.0;
		U[k + 1][0] = (d0 - c0 * U[k + 1][1]) / b0;
		U[k + 1][N - 1] = (dN - aN * U[k + 1][N - 2]) / bN;
	}
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

void implicit2Points1Order(vector<vector<double>>& U,
	vector<double>& lower, vector<double>& main, vector<double>& upper,
	vector<double>& coeffs, bool getA, int k) {
	double aj = a * tau / (h * h) - tau * b / 2.0 / h;
	double bj = tau * (c - 2.0 * a / (h * h)) - 1.0;
	double cj = a * tau / (h * h) + tau * b / 2.0 / h;
	if (getA) {
		double b0 = 1.0 - 1.0 / h;
		double c0 = 1.0 / h;
		double aN = -1.0 / h;
		double bN = 1.0 + 1.0 / h;
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
			coeffs[j] = -U[k][j];
		}
	}
}
void implicit3Points2Order(vector<vector<double>>& U,
	vector<double>& lower, vector<double>& main, vector<double>& upper,
	vector<double>& coeffs, bool getA, int k) {
	double aj = a * tau / (h * h) - tau * b / 2.0 / h;
	double bj = tau * (c - 2.0 * a / (h * h)) - 1.0;
	double cj = a * tau / (h * h) + tau * b / 2.0 / h;
	if (getA) {
		double a0 = 1.0 - 3.0 * 1.0 / h / 2.0;
		double b0 = 2.0 * 1.0 / h;
		double c0 = -1.0 / h / 2.0;
		double aN = 1.0 / h / 2.0;
		double bN = -2.0 * 1.0 / h;
		double cN = 1.0 + 3.0 * 1.0 / h / 2.0;

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
			coeffs[j] = -U[k][j];
		}
	}
}
void implicit2Points2Order(vector<vector<double>>& U,
	vector<double>& lower, vector<double>& main, vector<double>& upper,
	vector<double>& coeffs, bool getA, int k) {
	double aj = a * tau / (h * h) - tau * b / 2.0 / h;
	double bj = tau * (c - 2.0 * a / (h * h)) - 1.0;
	double cj = a * tau / (h * h) + tau * b / 2.0 / h;
	if (getA) {
		double b0 = 2.0 * a / h + h / tau - h * c - 1.0 / 1.0 * (2.0 * a - b * h);
		double c0 = -2.0 * a / h;
		double aN = -2.0 * a / h;
		double bN = 2.0 * a / h + h / tau - h * c + 1.0 / 1.0 * (2.0 * a + b * h);

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
	}
	else {
		coeffs[0] = h / tau * U[k][0] - phi0(tau * (k + 1.0)) * (2.0 * a - b * h) / 1.0;
		coeffs[N - 1] = h / tau * U[k][N - 1] + phi1(tau * (k + 1.0)) * (2.0 * a + b * h) / 1.0;
		for (int j = 1; j < N - 1; ++j) {
			coeffs[j] = -U[k][j];
		}
	}
}

vector<vector<double>> implicit2p1o() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = sin(j * h);
	}
	vector<double> lower(N);
	vector<double> main(N);
	vector<double> upper(N);
	vector<double> coeffs(N);
	implicit2Points1Order(U, lower, main, upper, coeffs, true, 0);
	for (int k = 0; k < K - 1; ++k) {
		implicit2Points1Order(U, lower, main, upper, coeffs, false, k);
		tridiagonalAlgo(lower, main, upper, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}

vector<vector<double>> implicit3p2o() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = sin(j * h);
	}
	vector<double> lower(N), main(N), upper(N);
	vector<double> coeffs(N);
	implicit3Points2Order(U, lower, main, upper, coeffs, true, 0);
	for (int k = 0; k < K - 1; ++k) {
		implicit3Points2Order(U, lower, main, upper, coeffs, false, k);
		tridiagonalAlgo(lower, main, upper, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}

vector<vector<double>> implicit2p2o() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = sin(j * h);
	}
	vector<double> lower(N);
	vector<double> main(N);
	vector<double> upper(N);
	vector<double> coeffs(N);
	implicit2Points2Order(U, lower, main, upper, coeffs, true, 0);
	for (int k = 0; k < K - 1; ++k) {
		implicit2Points2Order(U, lower, main, upper, coeffs, false, k);
		tridiagonalAlgo(lower, main, upper, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}

vector<vector<double>> crankNicolson() {
	vector<vector<double>> U(K, vector<double>(N, 0.0));
	vector<double> l(N), m(N), u(N), coeffs(N);
	for (int j = 0; j < N; ++j) {
		U[0][j] = sin(j * h);
	}
	double aj = sigma * a * tau / (h * h) - tau * b / 2.0 / h;
	double bj = tau * (c - sigma * 2.0 * a / (h * h)) - 1.0;
	double cj = sigma * a * tau / (h * h) + tau * b / 2.0 / h;
	double a0 = 1.0 - 3.0 * 1.0 / h / 2.0;
	double b0 = 2.0 * 1.0 / h;
	double c0 = -1.0 / h / 2.0;
	double aN = 1.0 / h / 2.0;
	double bN = -2.0 * 1.0 / h;
	double cN = 1.0 + 3.0 * 1.0 / h / 2.0;
	l[0] = a0;
	m[0] = b0;
	u[0] = c0;
	for (int j = 1; j < N - 1; ++j) {
		l[j] = aj;
		m[j] = bj;
		u[j] = cj;
	}
	l[N - 1] = aN;
	m[N - 1] = bN;
	u[N - 1] = cN;
	for (int k = 0; k < K - 1; ++k) {
		coeffs[0] = phi0(tau * (k + 1));
		coeffs[N - 1] = phi1(tau * (k + 1));
		for (int j = 1; j < N - 1; ++j) {
			coeffs[j] = -U[k][j] - (1.0 - sigma) * (a * tau / (h * h)) * (U[k][j + 1] - 2.0 * U[k][j] + U[k][j - 1]);
		}
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
		file << a << " " << b << " " << c << " " << sigma << " " << N << " " << timeSlice << endl;
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
	a = 12;
	b = 0.4;
	c = -19;
	sigma = 0.5;
	N = 100;
	timeSlice = 200;
	K = 1000;
	h = pi / N;
	tau = sigma * h * h / a;
	mu = b * tau / 2 / h;
	vector<vector<double>> U1, U2, U3, U4, U5, U6, U7;
	vector<string> fileNames{ "explicit2point1order.txt", 
                                "explicit3point2order.txt", 
                                "explicit2point2order.txt", 
                                "implicit2point1order.txt", 
                                "implicit3point2order.txt", 
                                "implicit2point2order.txt", 
                                "crankNicholson.txt" };
	U1 = explicit2Points1Order();
	U2 = explicit3Points2Order();
	U3 = explicit2Points2Order();
	U4 = implicit2p1o();
	U5 = implicit3p2o();
	U6 = implicit2p2o();
    U7 = crankNicolson();
    writeVectorToFile(U1, fileNames[0]);
	writeVectorToFile(U2, fileNames[1]);
	writeVectorToFile(U3, fileNames[2]);
	writeVectorToFile(U4, fileNames[3]);
	writeVectorToFile(U5, fileNames[4]);
	writeVectorToFile(U6, fileNames[5]);
	writeVectorToFile(U7, fileNames[6]);
}