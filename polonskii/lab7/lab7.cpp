/*
Polonskii Kirill, M8O-408B-20. Task 10 (20th in the group list)
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <numbers>
#include <string>

// Constants of task 10
const double L = 1.57079632679;
// Boundary conditions
double phi0(double y); // u(0, y)
double phi1(double y); // u(L, y)
double phi2(double x); // u(x, 0)
double phi3(double x); // u(x, L)
// Analytic solution
double analSol(double x, double y);
// Solution methods
std::vector<std::vector<double>> solveEllipEquation(int meth);
void interpolate(std::vector<std::vector<double>>& U);
void methLiebmann(std::vector<std::vector<double>>& U);
void methSeidel(std::vector<std::vector<double>>& U);
void methUpperRelax(std::vector<std::vector<double>>& U);
// Error function (comparison between U and analytical solution)
std::vector<double> getError(std::vector<std::vector<double>>& U);
// Inaccuracy function (comparison between old U and new U)
double getMaxInaccuracy(std::vector<std::vector<double>>& oldU, std::vector<std::vector<double>>& newU);

double a, b, c;
double omega, eps, Ny, sliceY, Nx;
double hX, hY;

int main() {

	std::cout << "Do you want to enter custom parameters? (y/n)\n";
	char enterInitParams;
	std::cin >> enterInitParams;

	if (enterInitParams == 'y') {
		std::cout << "Enter coefficients a, b, c separated by a space\n";
		std::cin >> a >> b >> c;
		std::cout << "Enter omega, eps (for methods) and sliceY separated by a space\n";
		std::cin >> omega >> eps >> sliceY;
	}
	else {
		a = -2.0;
		b = -2.0;
		c = -4.0;
		omega = 1.5;
		eps = 0.0001;
		sliceY = 10;
	}
	Ny = 100;
	Nx = 100;
	hX = L / Nx;
	hY = L / Ny;

	std::cout << "Choose a method:\n\
1. Liebmann method\n\
2. Seidel method\n\
3. Simple iterations with upper-relaxation\n\
";
	int method;
	std::cin >> method;
	std::vector<std::vector<double>> U;
	std::string fileName;
	U = solveEllipEquation(method);
	if (method == 1) {
		fileName = "Lieb.txt";
	}
	else if (method == 2) {
		fileName = "Seid.txt";
	}
	else if (method == 3) {
		fileName = "UpRe.txt";
	}
	
	std::cout << std::endl;
	std::vector<double> error = getError(U);
	std::cout << "Result will be written in file: " << fileName << std::endl;
	std::ofstream file(fileName);
	file << Ny << " " << Nx << " " << sliceY << std::endl;
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			file << U[i][j] << " ";
		}
		file << std::endl;
	}
	for (int i = 0; i < Nx; ++i) {
		file << error[i] << ' ';
	}
	file.close();

	return 0;
}

// Boundary conditions
double phi0(double y) { // u(0, y)
	return std::exp(-y) * std::cos(y);
}
double phi1(double y) { // u(L, y)
	return 0.0;
}
double phi2(double x) { // u(x, 0)
	return std::exp(-x) * std::cos(x);
}
double phi3(double x) { // u(x, L)
	return 0.0;
}
// Analytic solution
double analSol(double x, double y) {
	return std::exp(-x - y) * std::cos(x) * std::cos(y);
}
// Solution methods
std::vector<std::vector<double>> solveEllipEquation(int meth) {
	std::vector<std::vector<double>> U(Nx, std::vector<double>(Ny, 0.0));
	for (int i = 0; i < Nx; ++i) {
		U[i][0] = phi2(hX * i);
		U[i][Ny - 1] = phi3(hX * i);
	}
	for (int j = 0; j < Ny; ++j) {
		U[0][j] = phi0(hY * j);
		U[Nx - 1][j] = phi1(hY * j);
	}
	interpolate(U);
	if (meth == 1) {
		methLiebmann(U);
	}
	else if (meth == 2) {
		methSeidel(U);
	}
	else if (meth == 3) {
		methUpperRelax(U);
	}

	return U;
			
}
void interpolate(std::vector<std::vector<double>>& U) {
	for (int i = 1; i < Nx - 1; ++i) {
		for (int j = 1; j < Ny - 1; ++j) {
			double alpha = (j * hY) / L;
			U[i][j] = phi2(i * hX) * (1 - alpha) + phi3(i * hX) * alpha;
		}
	}
			
}
void methLiebmann(std::vector<std::vector<double>>& U) {
	double delta = 1 / (2 / (hX * hX) + 2 / (hY * hY) + c);
	double hXCoeff = 1 / (hX * hX);
	double aCoeff = a / 2 / hX;
	double hYCoeff = 1 / (hY * hY);
	double bCoeff = b / 2 / hY;
	int n = 0;
	while (true) {
		++n;
		std::vector<std::vector<double>> prevU = U;
		for (int i = 1; i < Nx - 1; ++i) {
			for (int j = 1; j < Ny - 1; ++j) {
				U[i][j] = delta * ((hXCoeff + aCoeff) * prevU[i - 1][j] +
					(hXCoeff - aCoeff) * prevU[i + 1][j] +
					(hYCoeff + bCoeff) * prevU[i][j - 1] +
					(hYCoeff - bCoeff) * prevU[i][j + 1]);
					
			}
		}
		double inaccuracy = getMaxInaccuracy(prevU, U);
		if (inaccuracy < eps)
			break;
	}
}
void methSeidel(std::vector<std::vector<double>>& U) {
	double delta = 1 / (2 / (hX * hX) + 2 / (hY * hY) + c);
	double hXCoeff = 1 / (hX * hX);
	double aCoeff = a / 2 / hX;
	double hYCoeff = 1 / (hY * hY);
	double bCoeff = b / 2 / hY;
	int n = 0;
	while (true) {
		++n;
		std::vector<std::vector<double>> prevU = U;
		for (int i = 1; i < Nx - 1; ++i) {
			for (int j = 1; j < Ny - 1; ++j) {
				U[i][j] = delta * ((hXCoeff + aCoeff) * U[i - 1][j] + (hXCoeff - aCoeff) * U[i + 1][j] +
					(hYCoeff + bCoeff) * U[i][j - 1] +
					(hYCoeff - bCoeff) * U[i][j + 1]);

			}
		}
		double inaccuracy = getMaxInaccuracy(prevU, U);
		if (inaccuracy < eps)
			break;
	}
}
void methUpperRelax(std::vector<std::vector<double>>& U) {
	double delta = 1 / (2 / (hX * hX) + 2 / (hY * hY) + c);
	double hXCoeff = 1 / (hX * hX);
	double aCoeff = a / 2 / hX;
	double hYCoeff = 1 / (hY * hY);
	double bCoeff = b / 2 / hY;
	int n = 0;
	bool converge = false;
	while (n < 1000) {
		++n;
		std::vector<std::vector<double>> prevU = U;
		for (int i = 1; i < Nx - 1; ++i) {
			for (int j = 1; j < Ny - 1; ++j) {
				U[i][j] = U[i][j] + omega * (delta * ((hXCoeff + aCoeff) * U[i - 1][j] +
					(hXCoeff - aCoeff) * prevU[i + 1][j] +
					(hYCoeff + bCoeff) * U[i][j - 1] +
					(hYCoeff - bCoeff) * prevU[i][j + 1]) - U[i][j]);
			}
		}
		double inaccuracy = getMaxInaccuracy(prevU, U);
		if (inaccuracy < eps) {
			converge = true;
			break;
		}
			
	}
	if (!converge) {
		std::cout << "Method don't converge\n";
	}
}
// Error function
std::vector<double> getError(std::vector<std::vector<double>>& U) {
	std::vector<double> error(Nx, 0.0);
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			double curError = std::abs(U[i][j] - analSol(i * hX, j * hY));
			if (curError > error[i]) {
				error[i] = curError;
			}
		}
	}
	return error;
}
// Inaccuracy function (comparison between old U and new U)
double getMaxInaccuracy(std::vector<std::vector<double>>& oldU, std::vector<std::vector<double>>& newU) {
	double maxInaccuracy = -1.0;
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			double curInaccuracy = std::abs(newU[i][j] - oldU[i][j]);
			if (curInaccuracy > maxInaccuracy) {
				maxInaccuracy = curInaccuracy;
			}
		}
	}
	return maxInaccuracy;
}
