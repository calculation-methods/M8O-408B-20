/*
Polonskii Kirill, M8O-408B-20. Task 10 (20th in the group list)
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <numbers>
#include <string>

// Constants of task 10
const double ALPHA0 = 1.0; // coeff before du/dx(0, t)
const double BETA0 = 1.0; // coeff before u(0, t)
const double ALPHA1 = 1.0; // coeff before du/dx(pi, t)
const double BETA1 = 1.0; // coeff before u(pi, t)
const double L = 3.14159265359;

// Boundary conditions
double phi0(double t);
double phi1(double t);
// Initial condition
double psi(double x);
// Analytic solution
double analSol(double x, double t);
// Explicit finite difference
std::vector<std::vector<double>> explicitMethod(int approx);
void explicit2Points1Order(std::vector<std::vector<double>>& U);
void explicit3Points2Order(std::vector<std::vector<double>>& U);
void explicit2Points2Order(std::vector<std::vector<double>>& U);
// Thomas algorithm
double tridiagonalAlgo(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d,
	std::vector<double>& x, int step, double prevP, double prevQ);
// Implicit finite difference
std::vector<std::vector<double>> implicitMethod(int approx);
void implicit2Points1Order(std::vector<std::vector<double>>& U, 
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k);
void implicit3Points2Order(std::vector<std::vector<double>>& U,
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k);
void implicit2Points2Order(std::vector<std::vector<double>>& U,
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k);
// Crank-Nicolson methodа
std::vector<std::vector<double>> CrankNicolsonMethod();
// Error function
std::vector<double> getError(std::vector<std::vector<double>>& U);

double a, b, c;
double sigma, N, timeSlice, K;
double h, tau, mu;

int main() {
	std::cout << "Do you want to enter custom parameters? (y/n)\n";
	char enterInitParams;
	std::cin >> enterInitParams;
	
	if (enterInitParams == 'y') {
		std::cout << "Enter coefficients a>0, b>0, c<0 separated by a space\n";
		std::cin >> a >> b >> c;
		std::cout << "Enter grid parameters sigma, N, timeSlice separated by a space\n";
		std::cin >> sigma >> N >> timeSlice;
	}
	else {
		a = 12;
		b = 0.4;
		c = -19;
		sigma = 0.5;
		N = 100;
		timeSlice = 200;
	}
	K = 1000;
	h = L / N;
	tau = sigma * h * h / a;
	mu = b * tau / 2 / h;

	std::cout << "Choose a method:\n\
1. Explicit finite difference with 2 points 1 order approximation\n\
2. Explicit finite difference with 3 points 2 order approximation\n\
3. Explicit finite difference with 2 points 2 order approximation\n\
4. Implicit finite difference with 2 points 1 order approximation\n\
5. Implicit finite difference with 3 points 2 order approximation\n\
6. Implicit finite difference with 2 points 2 order approximation\n\
7. Crank-Nicolson scheme\n\
";
	int method;
	std::cin >> method;
	std::vector<std::vector<double>> U;
	std::string fileName;
	if (method == 1) {
		U = explicitMethod(1);
		fileName = "EFD_2P1O.txt";
	}
	else if (method == 2) {
		U = explicitMethod(2);
		fileName = "EFD_3P2O.txt";
	}
	else if (method == 3) {
		U = explicitMethod(3);
		fileName = "EFD_2P2O.txt";
	}
	else if (method == 4) {
		U = implicitMethod(1);
		fileName = "IFD_2P1O.txt";
	}
	else if (method == 5) {
		U = implicitMethod(2);
		fileName = "IFD_3P2O.txt";
	}
	else if (method == 6) {
		U = implicitMethod(3);
		fileName = "IFD_2P2O.txt";
	}
	else if (method == 7) {
		U = CrankNicolsonMethod();
		fileName = "CNS.txt";
	}
	std::cout << std::endl;
	std::vector<double> error = getError(U);
	std::cout << "Result will be written in file: " << fileName << std::endl;
	std::ofstream file(fileName);
	file << a << " " << b << " " << c << " " << sigma << " " << N << " " << timeSlice << std::endl;
	for (int k = 0; k < K; ++k) {
		for (int j = 0; j < N; ++j) {
			file << U[k][j] << " ";
		}
		file << std::endl;
	}
	for (int k = 0; k < K; ++k) {
		file << error[k] << ' ';
	}
	file.close();
}

// Boundary conditions
double phi0(double t) {
	return std::exp((c - a) * t) * (std::cos(b * t) + std::sin(b * t));
}
double phi1(double t) {
	return -phi0(t);
}
// Initial condition
double psi(double x) {
	return std::sin(x);
}
// Analytic solution
double analSol(double x, double t) {
	return std::exp((c - a) * t) * std::sin(x + b * t);
}

// Explicit finite difference
std::vector<std::vector<double>> explicitMethod(int approx) {
	std::vector<std::vector<double>> U(K, std::vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi(j * h);
	}
	for (int k = 0; k < K - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			U[k + 1][j] = (sigma + mu) * U[k][j + 1] + (1 - 2 * sigma + tau * c) * U[k][j] + (sigma - mu) * U[k][j - 1];
		}
	}
	if (approx == 1) {
		explicit2Points1Order(U);
	}
	else if(approx == 2) {
		explicit3Points2Order(U);
	}
	else if (approx == 3) {
		explicit2Points2Order(U);
	}
	return U;
}
void explicit2Points1Order(std::vector<std::vector<double>>& U) {
	for (int k = 0; k < K - 1; ++k) {
		U[k + 1][0] = (-ALPHA0 / h * U[k + 1][1] + phi0(tau * (k + 1))) / (BETA0 - ALPHA0 / h);
		U[k + 1][N - 1] = (ALPHA1 / h * U[k + 1][N - 2] + phi1(tau * (k + 1))) / (BETA1 + ALPHA1 / h);
	}
}
void explicit3Points2Order(std::vector<std::vector<double>>& U) {
	for (int k = 0; k < K - 1; ++k) {
		U[k + 1][0] = (-ALPHA0 / h / 2 * (4 * U[k + 1][1] - U[k + 1][2]) + phi0(tau * (k + 1))) /
					(BETA0 - 3 * ALPHA0 / h / 2);
		U[k + 1][N - 1] = (-ALPHA1 / h / 2 * (U[k + 1][N - 3] - 4 * U[k + 1][N - 2]) + phi1(tau * (k + 1))) /
						(BETA1 + 3 * ALPHA1 / h / 2);
	}
}
void explicit2Points2Order(std::vector<std::vector<double>>& U) {
	double b0 = 2.0 * a / h + h / tau - h * c - BETA0 / ALPHA0 * (2 * a - b * h);
	double c0 = -2.0 * a / h;
	double bN = 2.0 * a / h + h / tau - h * c + BETA1 / ALPHA1 * (2 * a + b * h);
	double aN = -2.0 * a / h;
	for (int k = 0; k < K - 1; ++k) {
		double d0 = h / tau * U[k][0] - phi0(tau * (k + 1)) * (2 * a - b * h) / ALPHA0;
		double dN = h / tau * U[k][N - 1] + phi1(tau * (k + 1)) * (2 * a + b * h) / ALPHA1;
		U[k + 1][0] = (d0 - c0 * U[k + 1][1]) / b0;
		U[k + 1][N - 1] = (dN - aN * U[k + 1][N - 2]) / bN;
	}
}

// Thomas algorithm
double tridiagonalAlgo(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d,
	std::vector<double>& x, int step, double prevP, double prevQ) {
	if (step == a.size()) {
		return 0;
	}
	double p = ((-1.0) * c[step]) / (b[step] + a[step] * prevP);
	double q = (d[step] - a[step] * prevQ) / (b[step] + a[step] * prevP);
	x[step] = p * tridiagonalAlgo(a, b, c, d, x, step + 1, p, q) + q;
	return x[step];
}

// Implicit finite difference
std::vector<std::vector<double>> implicitMethod(int approx) {
	std::vector<std::vector<double>> U(K, std::vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi(j * h);
	}
	// get only once 3 diagonals in matrix A
	std::vector<double> lower(N);
	std::vector<double> main(N);
	std::vector<double> upper(N);
	std::vector<double> coeffs(N);
	if (approx == 1) {
		implicit2Points1Order(U, lower, main, upper, coeffs, true, 0);
	}
	else if (approx == 2) {
		implicit3Points2Order(U, lower, main, upper, coeffs, true, 0);
	}
	else {
		implicit2Points2Order(U, lower, main, upper, coeffs, true, 0);
	}
	for (int k = 0; k < K - 1; ++k) {
		if (approx == 1) {
			implicit2Points1Order(U, lower, main, upper, coeffs, false, k);
		}
		else if (approx == 2) {
			implicit3Points2Order(U, lower, main, upper, coeffs, false, k);
		}
		else {
			implicit2Points2Order(U, lower, main, upper, coeffs, false, k);
		}
		tridiagonalAlgo(lower, main, upper, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}
void implicit2Points1Order(std::vector<std::vector<double>>& U,
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k) {
	double aj = a * tau / (h * h) - tau * b / 2.0 / h;
	double bj = tau * (c - 2.0 * a / (h * h)) - 1.0;
	double cj = a * tau / (h * h) + tau * b / 2.0 / h;
	if (getA) { // get only once 3 diagonals in matrix A
		double b0 = BETA0 - ALPHA0 / h;
		double c0 = ALPHA0 / h;
		double aN = -ALPHA1 / h;
		double bN = BETA1 + ALPHA1 / h;
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
	else { // get system coefficients
		coeffs[0] = phi0(tau * (k + 1));
		coeffs[N - 1] = phi1(tau * (k + 1));
		for (int j = 1; j < N - 1; ++j) {
			coeffs[j] = -U[k][j];
		}
	}
}
void implicit3Points2Order(std::vector<std::vector<double>>& U,
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k) {
	double aj = a * tau / (h * h) - tau * b / 2.0 / h;
	double bj = tau * (c - 2.0 * a / (h * h)) - 1.0;
	double cj = a * tau / (h * h) + tau * b / 2.0 / h;
	if (getA) { // get only once 3 diagonals in matrix A
		double a0 = BETA0 - 3.0 * ALPHA0 / h / 2.0;
		double b0 = 2.0 * ALPHA0 / h;
		double c0 = -ALPHA0 / h / 2.0;
		double aN = ALPHA1 / h / 2.0;
		double bN = -2.0 * ALPHA1 / h;
		double cN = BETA1 + 3.0 * ALPHA1 / h / 2.0;

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
	else { // get system coefficients
		coeffs[0] = phi0(tau * (k + 1));
		coeffs[N - 1] = phi1(tau * (k + 1));
		for (int j = 1; j < N - 1; ++j) {
			coeffs[j] = -U[k][j];
		}
	}
}
void implicit2Points2Order(std::vector<std::vector<double>>& U,
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k) {
	double aj = a * tau / (h * h) - tau * b / 2.0 / h;
	double bj = tau * (c - 2.0 * a / (h * h)) - 1.0;
	double cj = a * tau / (h * h) + tau * b / 2.0 / h;
	if (getA) { // get only once 3 diagonals in matrix A
		double b0 = 2.0 * a / h + h / tau - h * c - BETA0 / ALPHA0 * (2.0 * a - b * h);
		double c0 = -2.0 * a / h;
		double aN = -2.0 * a / h;
		double bN = 2.0 * a / h + h / tau - h * c + BETA1 / ALPHA1 * (2.0 * a + b * h);

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
	else { // get system coefficients
		coeffs[0] = h / tau * U[k][0] - phi0(tau * (k + 1.0)) * (2.0 * a - b * h) / ALPHA0;
		coeffs[N - 1] = h / tau * U[k][N - 1] + phi1(tau * (k + 1.0)) * (2.0 * a + b * h) / ALPHA1;
		for (int j = 1; j < N - 1; ++j) {
			coeffs[j] = -U[k][j];
		}
	}
}

// Crank-Nicolson methodа
std::vector<std::vector<double>> CrankNicolsonMethod() {
	std::vector<std::vector<double>> U(K, std::vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi(j * h);
	}
	double aj = sigma * a * tau / (h * h) - tau * b / 2.0 / h;
	double bj = tau * (c - sigma * 2.0 * a / (h * h)) - 1.0;
	double cj = sigma * a * tau / (h * h) + tau * b / 2.0 / h;
	// get only once 3 diagonals in matrix A
	std::vector<double> lower(N);
	std::vector<double> main(N);
	std::vector<double> upper(N);
	std::vector<double> coeffs(N);

	double a0 = BETA0 - 3.0 * ALPHA0 / h / 2.0;
	double b0 = 2.0 * ALPHA0 / h;
	double c0 = -ALPHA0 / h / 2.0;
	double aN = ALPHA1 / h / 2.0;
	double bN = -2.0 * ALPHA1 / h;
	double cN = BETA1 + 3.0 * ALPHA1 / h / 2.0;

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
	for (int k = 0; k < K - 1; ++k) {
		coeffs[0] = phi0(tau * (k + 1));
		coeffs[N - 1] = phi1(tau * (k + 1));
		for (int j = 1; j < N - 1; ++j) {
			coeffs[j] = -U[k][j] - (1.0 - sigma) * (a * tau / (h * h)) * (U[k][j + 1] - 2.0 * U[k][j] + U[k][j - 1]);
		}
		tridiagonalAlgo(lower, main, upper, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}

// Error function
std::vector<double> getError(std::vector<std::vector<double>>& U) {
	std::vector<double> error(K, 0.0);
	for (int k = 0; k < K; ++k) {
		for (int j = 0; j < N; ++j) {
			double curError = std::abs(U[k][j] - analSol(h * j, tau * k));
			if (curError > error[k]) {
				error[k] = curError;
			}
		}
	}
	return error;
}
