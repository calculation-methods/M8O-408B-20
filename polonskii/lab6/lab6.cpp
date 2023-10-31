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
const double BETA0 = 0.0; // coeff before u(0, t)
const double ALPHA1 = 1.0; // coeff before du/dx(pi, t)
const double BETA1 = 0.0; // coeff before u(pi, t)
const double L = 3.14159265359;

// Boundary conditions
double phi0(double t);
double phi1(double t);
// Initial conditions
double psi0(double x);
double psi1(double x);
// Function in task
double func(double x, double t);
// Analytic solution
double analSol(double x, double t);
// Explicit finite difference
std::vector<std::vector<double>> explicitMethod(int initApprox, int boundApprox); // changed
void explicitBound2Points1Order(std::vector<std::vector<double>>& U); // changed
void explicitBound3Points2Order(std::vector<std::vector<double>>& U); // changed
void explicitBound2Points2Order(std::vector<std::vector<double>>& U); // changed
void explicitInit1Order(std::vector<std::vector<double>>& U); // changed
void explicitInit2Order(std::vector<std::vector<double>>& U); // changed
// Thomas algorithm
double tridiagonalAlgo(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d,
	std::vector<double>& x, int step, double prevP, double prevQ);
// Implicit finite difference
std::vector<std::vector<double>> implicitMethod(int initApprox, int boundApprox); // changed
void implicitBound2Points1Order(std::vector<std::vector<double>>& U,
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k); // changed
void implicitBound3Points2Order(std::vector<std::vector<double>>& U,
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k); // changed
void implicitBound2Points2Order(std::vector<std::vector<double>>& U,
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k); // changed
void implicitInit1Order(std::vector<std::vector<double>>& U); // changed
void implicitInit2Order(std::vector<std::vector<double>>& U); // changed
// Error function
std::vector<double> getError(std::vector<std::vector<double>>& U);
// Derivatives
double deriv(double (*f)(double), double x);
double deriv2(double (*innerFunc)(double x1), double x);
	

double a, b, c, d;
double sigma, N, timeSlice, K;
double h, tau, mu;

int main() {
	std::cout << "Do you want to enter custom parameters? (y/n)\n";
	char enterInitParams;
	std::cin >> enterInitParams;

	if (enterInitParams == 'y') {
		std::cout << "Enter coefficients a>0, b>0, c<0, d>0 separated by a space\n";
		std::cin >> a >> b >> c >> d;
		std::cout << "Enter grid parameters sigma, N, timeSlice separated by a space\n";
		std::cin >> sigma >> N >> timeSlice;
	}
	else {
		a = 1.0;
		b = 1.0;
		c = -1.0;
		d = 3.0;
		sigma = 0.5;
		N = 100;
		timeSlice = 20;
	}
	K = 500;
	h = L / N;
	tau = std::sqrt(sigma * h * h / a);
	mu = b * tau * tau / 2.0 / h;

	std::cout << "Choose a method:\n\
1. Explicit finite difference: 2 points 1 order approx for boundary conditions, 1 order approx for initial conditions\n\
2. Explicit finite difference: 2 points 1 order approx for boundary conditions, 2 order approx for initial conditions\n\
3. Explicit finite difference: 3 points 2 order approx for boundary conditions, 1 order approx for initial conditions\n\
4. Explicit finite difference: 3 points 2 order approx for boundary conditions, 2 order approx for initial conditions\n\
5. Explicit finite difference: 2 points 2 order approx for boundary conditions, 1 order approx for initial conditions\n\
6. Explicit finite difference: 2 points 2 order approx for boundary conditions, 2 order approx for initial conditions\n\
7. Implicit finite difference: 2 points 1 order approx for boundary conditions, 1 order approx for initial conditions\n\
8. Implicit finite difference: 2 points 1 order approx for boundary conditions, 2 order approx for initial conditions\n\
9. Implicit finite difference: 3 points 2 order approx for boundary conditions, 1 order approx for initial conditions\n\
10. Implicit finite difference: 3 points 2 order approx for boundary conditions, 2 order approx for initial conditions\n\
11. Implicit finite difference: 2 points 2 order approx for boundary conditions, 1 order approx for initial conditions\n\
12. Implicit finite difference: 2 points 2 order approx for boundary conditions, 2 order approx for initial conditions\n\
";
	int method;
	std::cin >> method;
	std::vector<std::vector<double>> U;
	std::string fileName;
	if (method == 1) {
		U = explicitMethod(1, 1);
		fileName = "EFD_1O_2P1O.txt";
	}
	else if (method == 2) {
		U = explicitMethod(2, 1);
		fileName = "EFD_2O_2P1O.txt";
	}
	else if (method == 3) {
		U = explicitMethod(1, 2);
		fileName = "EFD_1O_3P2O.txt";
	}
	else if (method == 4) {
		U = explicitMethod(2, 2);
		fileName = "EFD_2O_3P2O.txt";
	}
	else if (method == 5) {
		U = explicitMethod(1, 3);
		fileName = "EFD_1O_2P2O.txt";
	}
	else if (method == 6) {
		U = explicitMethod(2, 3);
		fileName = "EFD_2O_2P2O.txt";
	}
	else if (method == 7) {
		U = implicitMethod(1, 1);
		fileName = "IFD_1O_2P1O.txt";
	}
	else if (method == 8) {
		U = implicitMethod(2, 1);
		fileName = "IFD_2O_2P1O.txt";
	}
	else if (method == 9) {
		U = implicitMethod(1, 2);
		fileName = "IFD_1O_3P2O.txt";
	}
	else if (method == 10) {
		U = implicitMethod(2, 2);
		fileName = "IFD_2O_3P2O.txt";
	}
	else if (method == 11) {
		U = implicitMethod(1, 3);
		fileName = "IFD_1O_2P2O.txt";
	}
	else if (method == 12) {
		U = implicitMethod(2, 3);
		fileName = "IFD_2O_2P2O.txt";
	}
	std::cout << std::endl;
	std::vector<double> error = getError(U);
	std::cout << "Result will be written in file: " << fileName << std::endl;
	std::ofstream file(fileName);
	file << a << " " << b << " " << c << " " << d << " " << sigma << " " << N << " " << timeSlice << std::endl;
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
	return std::exp(-t);
}
double phi1(double t) {
	return -phi0(t);
}
// Initial conditions
double psi0(double x) {
	return std::sin(x);
}
double psi1(double x) {
	return -std::sin(x);
}
// Function in task
double func(double x, double t) {
	return -std::cos(x) * std::exp(-t);
}
// Analytic solution
double analSol(double x, double t) {
	return std::exp(-t) * std::sin(x);
}

// Explicit finite difference
std::vector<std::vector<double>> explicitMethod(int initApprox, int boundApprox) {
	std::vector<std::vector<double>> U(K, std::vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
	if (initApprox == 1) {
		explicitInit1Order(U);
	}
	else if (initApprox == 2) {
		explicitInit1Order(U);
	}
	
	for (int k = 1; k < K - 1; ++k) {
		for (int j = 1; j < N - 1; ++j) {
			U[k + 1][j] = ((sigma + mu) * U[k][j + 1] + (-2 * sigma + 2 + c * tau * tau) * U[k][j] + (sigma - mu) * U[k][j - 1]
				+ (-1 + d * tau / 2) * U[k - 1][j] + tau * tau * func(j * h, k * tau)) / (1 + d * tau / 2);
		}
	}
	if (boundApprox == 1) {
		explicitBound2Points1Order(U);
	}
	else if (boundApprox == 2) {
		explicitBound3Points2Order(U);
	}
	else if (boundApprox == 3) {
		explicitBound2Points2Order(U);
	}
	return U;
}
void explicitBound2Points1Order(std::vector<std::vector<double>>& U) {
	for (int k = 1; k < K - 1; ++k) {
		U[k + 1][0] = (-ALPHA0 / h * U[k + 1][1] + phi0(tau * (k + 1))) / (BETA0 - ALPHA0 / h);
		U[k + 1][N - 1] = (ALPHA1 / h * U[k + 1][N - 2] + phi1(tau * (k + 1))) / (BETA1 + ALPHA1 / h);
	}
}
void explicitBound3Points2Order(std::vector<std::vector<double>>& U) {
	for (int k = 1; k < K - 1; ++k) {
		U[k + 1][0] = (-ALPHA0 / h / 2 * (4 * U[k + 1][1] - U[k + 1][2]) + phi0(tau * (k + 1))) /
			(BETA0 - 3 * ALPHA0 / h / 2);
		U[k + 1][N - 1] = (-ALPHA1 / h / 2 * (U[k + 1][N - 3] - 4 * U[k + 1][N - 2]) + phi1(tau * (k + 1))) /
			(BETA1 + 3 * ALPHA1 / h / 2);
	}
}
void explicitBound2Points2Order(std::vector<std::vector<double>>& U) {
	double b0 = 2.0 * a / h + h / (tau * tau) - h * c - BETA0 / ALPHA0 * (2 * a - b * h) + d * h / 2 / tau;
	double c0 = -2.0 * a / h;
	double bN = 2.0 * a / h + h / (tau * tau) - h * c + BETA1 / ALPHA1 * (2 * a + b * h) + d * h / 2 / tau;
	double aN = -2.0 * a / h;
	for (int k = 1; k < K - 1; ++k) {
		double d0 = 2 * h / (tau * tau) * U[k][0] + (-h / (tau * tau) + d * h / 2 / tau) * U[k - 1][0] -
			phi0(tau * (k + 1)) * (2 * a - b * h) / ALPHA0 + h * func(0, tau * (k + 1));
		double dN = 2 * h / (tau * tau) * U[k][N - 1] + (-h / (tau * tau) + d * h / 2 / tau) * U[k - 1][N - 1] +
			phi1(tau * (k + 1)) * (2 * a + b * h) / ALPHA1 - h * func(L, tau * (k + 1));
		U[k + 1][0] = (d0 - c0 * U[k + 1][1]) / b0;
		U[k + 1][N - 1] = (dN - aN * U[k + 1][N - 2]) / bN;
	}
}
void explicitInit1Order(std::vector<std::vector<double>>& U) {
	for (int j = 0; j < N; ++j) {
		U[1][j] = psi0(j * h) + psi1(j * h) * tau;
	}
}
void explicitInit2Order(std::vector<std::vector<double>>& U) {
	for (int j = 0; j < N; ++j) {
		U[1][j] = psi0(j * h) + psi1(j * h) * (tau - d * (tau * tau) / 2) + (a *  deriv2(&psi0, j * h) + b * deriv(&psi0, j * h) +
			c * psi0(j * h) + func(j * h, tau)) * (tau * tau) / 2;
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
std::vector<std::vector<double>> implicitMethod(int initApprox, int boundApprox) {
	std::vector<std::vector<double>> U(K, std::vector<double>(N, 0.0));
	for (int j = 0; j < N; ++j) {
		U[0][j] = psi0(j * h);
	}
	if (initApprox == 1) {
		implicitInit1Order(U);
	}
	else {
		implicitInit2Order(U);
	}
	// get only once 3 diagonals in matrix A
	std::vector<double> lower(N);
	std::vector<double> main(N);
	std::vector<double> upper(N);
	std::vector<double> coeffs(N);
	if (boundApprox == 1) {
		implicitBound2Points1Order(U, lower, main, upper, coeffs, true, 0);
	}
	else if (boundApprox == 2) {
		implicitBound3Points2Order(U, lower, main, upper, coeffs, true, 0);
	}
	else {
		implicitBound2Points2Order(U, lower, main, upper, coeffs, true, 0);
	}
	for (int k = 1; k < K - 1; ++k) {
		if (boundApprox == 1) {
			implicitBound2Points1Order(U, lower, main, upper, coeffs, false, k);
		}
		else if (boundApprox == 2) {
			implicitBound3Points2Order(U, lower, main, upper, coeffs, false, k);
		}
		else {
			implicitBound2Points2Order(U, lower, main, upper, coeffs, false, k);
		}
		tridiagonalAlgo(lower, main, upper, coeffs, U[k + 1], 0.0, 0.0, 0.0);
	}
	return U;
}
void implicitBound2Points1Order(std::vector<std::vector<double>>& U,
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k) {
	double aCoeff = a / (h * h);
	double bCoeff = b / 2.0 / h;
	double tCoeff = 1.0 / (tau * tau);
	double dCoeff = d / 2.0 / tau;

	double aj = bCoeff - aCoeff;
	double bj = tCoeff + dCoeff + 2.0 * aCoeff;
	double cj = -bCoeff - aCoeff;
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
			coeffs[j] = U[k - 1][j] * (-tCoeff + dCoeff) + U[k][j] * (2.0 * tCoeff + c) + func(j * h, k * tau);
		}
	}
}
void implicitBound3Points2Order(std::vector<std::vector<double>>& U,
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k) {
	double aCoeff = a / (h * h);
	double bCoeff = b / 2.0 / h;
	double tCoeff = 1.0 / (tau * tau);
	double dCoeff = d / 2.0 / tau;

	double aj = bCoeff - aCoeff;
	double bj = tCoeff + dCoeff + 2.0 * aCoeff;
	double cj = -bCoeff - aCoeff;
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
			coeffs[j] = U[k - 1][j] * (-tCoeff + dCoeff) + U[k][j] * (2 * tCoeff + c) + func(j * h, k * tau);
		}
	}
}
void implicitBound2Points2Order(std::vector<std::vector<double>>& U,
	std::vector<double>& lower, std::vector<double>& main, std::vector<double>& upper,
	std::vector<double>& coeffs, bool getA, int k) {
	double aCoeff = a / (h * h);
	double bCoeff = b / 2.0 / h;
	double tCoeff = 1.0 / (tau * tau);
	double dCoeff = d / 2.0 / tau;

	double aj = bCoeff - aCoeff;
	double bj = tCoeff + dCoeff + 2.0 * aCoeff;
	double cj = -bCoeff - aCoeff;
	if (getA) { // get only once 3 diagonals in matrix A
		double b0 = 2.0 * a / h + h / (tau * tau) - h * c - BETA0 / ALPHA0 * (2.0 * a - b * h) + d * h / 2 / tau;
		double c0 = -2.0 * a / h;
		double aN = -2.0 * a / h;
		double bN = 2.0 * a / h + h / (tau * tau) - h * c + BETA1 / ALPHA1 * (2.0 * a + b * h) + d * h / 2 / tau;

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
		coeffs[0] = 2 * h / (tau * tau) * U[k][0] + (-h / (tau * tau) + d * h / 2 / tau) * U[k - 1][0] - phi0(tau * (k + 1)) * (
			2 * a - b * h) / ALPHA0 + h * func(0, tau * (k + 1));
		coeffs[N - 1] = 2 * h / (tau * tau) * U[k][N - 1] + (-h / (tau * tau) + d * h / 2 / tau) * U[k - 1][N - 1] + phi1(
			tau * (k + 1)) * (2 * a + b * h) / ALPHA1 - h * func(L, tau * (k + 1));
		for (int j = 1; j < N - 1; ++j) {
			coeffs[j] = U[k - 1][j] * (-tCoeff + dCoeff) + U[k][j] * (2 * tCoeff + c) + func(j * h, k * tau);
		}
	}
}
void implicitInit1Order(std::vector<std::vector<double>>& U) {
	for (int j = 0; j < N; ++j) {
		U[1][j] = psi0(j * h) + psi1(j * h) * tau;
	}
}
void implicitInit2Order(std::vector<std::vector<double>>& U) {
	for (int j = 0; j < N; ++j) {
		U[1][j] = psi0(j * h) + psi1(j * h) * (tau - d * (tau * tau) / 2) + (a * deriv2(&psi0, j * h) + b * deriv(&psi0, j * h) +
			c * psi0(j * h) + func(j * h, tau)) * (tau * tau) / 2;
	}
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
// Derivative
double deriv(double (*f)(double), double x) {
	double dx = 0.000001;
	return (f(x + dx) - f(x)) / dx;
}
double deriv2(double (*innerFunc)(double x1), double x) {
	double dx = 0.000001;
	return (deriv(innerFunc, x + dx) - deriv(innerFunc, x)) / dx;
}
