/*
Polonskii Kirill, M8O-408B-20. Task 1
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

typedef std::vector<std::vector<std::vector<double>>> gridF;
// Constants of task 1
const double PI = 3.14159265359;
// Boundary conditions
double phi0(double y, double t); // u(0, y, t)
double phi1(double y, double t); // u(PI, y, t)
double phi2(double x, double t); // u(x, 0, t)
double phi3(double x, double t); // u(x, PI, t)
// Initial condition
double psi(double x, double y); // u(x, y, 0)
// Analytic solution
double analSol(double x, double y, double t);
// Solution methods
// Thomas algorithm
double tridiagonalAlgo(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d,
	std::vector<double>& x, int step, double prevP, double prevQ);
// Alternating direction method
gridF alterDirec();
// Fractional steps method
gridF fractSteps();
// Error function
std::vector<double> getError(gridF& U);

double a, mu1, mu2;
double xRange[2] = { 0.0, PI };
double yRange[2] = { 0.0, PI };
double tRange[2] = { 0.0, 1.0 };
double hX = 0.01;
double hY = 0.01;
double tau = 0.01;
int Nt = (tRange[1] - tRange[0]) / tau;
int Nx = (xRange[1] - xRange[0]) / hX;
int Ny = (yRange[1] - yRange[0]) / hY;
int sliceT, sliceX;

int main() {

	std::cout << "Do you want to enter custom parameters? (y/n)\n";
	char enterInitParams;
	std::cin >> enterInitParams;

	if (enterInitParams == 'y') {
		std::cout << "Enter coefficients a, mu1, mu2 separated by a space\n";
		std::cin >> a >> mu1 >> mu2;
		std::cout << "Enter sliceT and sliceX separated by a space\n";
		std::cin >> sliceT >> sliceX;
	}
	else {
		a = 1.0;
		mu1 = 2.0;
		mu2 = 1.0;
		sliceT = 1;
		sliceX = 1;
	}

	std::cout << "Choose a method:\n\
1. Alternating direction method\n\
2. Fractional steps method\n\
";
	int method;
	std::cin >> method;
	gridF U;
	std::string fileName;
	
	if (method == 1) {
		fileName = "AlDi.txt";
		U = alterDirec();
	}
	else if (method == 2) {
		fileName = "FrSt.txt";
		U = fractSteps();
	}
	std::vector<double> error = getError(U);
	std::cout << "Result will be written in file: " << fileName << std::endl;
	std::ofstream file(fileName);
	file << a << " " << mu1 << " " << mu2 << " " << sliceT << " " << sliceX << std::endl;
	for (int y = 0; y < Ny; ++y) {
		file << U[sliceT][sliceX][y] << " ";
	}
	file << std::endl;
	for (int t = 0; t < Nt; ++t) {
		file << error[t] << ' ';
	}
	file.close();

	return 0;
}

// Boundary conditions
double phi0(double y, double t) { // u(0, y, t)
	return std::cos(mu2 * y) * std::exp(-(mu1 * mu1 + mu2 * mu2) * a * t);
}
double phi1(double y, double t) { // u(PI, y, t)
	return std::pow(-1, mu1) * std::cos(mu2 * y) * std::exp(-(mu1 * mu1 + mu2 * mu2) * a * t);
}
double phi2(double x, double t) { // u(x, 0, t)
	return std::cos(mu1 * x) * std::exp(-(mu1 * mu1 + mu2 * mu2) * a * t);
}
double phi3(double x, double t) { // u(x, PI, t)
	return std::pow(-1, mu2) * std::cos(mu1 * x) * std::exp(-(mu1 * mu1 + mu2 * mu2) * a * t);
}
// Initial condition
double psi(double x, double y) { // u(x, y, 0)
	return std::cos(mu1 * x) * std::cos(mu2 * y);
}
// Analytic solution
double analSol(double x, double y, double t) {
	return std::cos(mu1 * x) * std::cos(mu2 * y) * std::exp(-(mu1 * mu1 + mu2 * mu2) * a * t);
}
// Solution methods
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
gridF alterDirec() {
	gridF U = std::vector<std::vector<std::vector<double>>>(Nt,
		std::vector<std::vector<double>>(Nx, std::vector<double>(Ny, 0.0)));

	for (int x = 0; x < Nx; ++x) {
		for (int y = 0; y < Ny; ++y) {
			U[0][x][y] = psi(x * hX, y * hY);
		}
	}

	for (int t = 1; t < Nt; ++t) {
		std::vector<std::vector<double>> halftimeU = std::vector<std::vector<double>>(Nx, std::vector<double>(Ny, 0.0));

		for (int x = 0; x < Nx; ++x) {
			U[t][x][0] = phi2(x * hX, t * tau);
			U[t][x][Ny-1] = phi3(x * hX, t * tau);
			halftimeU[x][0] = phi2(x * hX, t * tau - tau / 2.0);
			halftimeU[x][Ny-1] = phi3(x * hX, t * tau - tau / 2.0);
		}
		for (int y = 0; y < Ny; ++y) {
			U[t][0][y] = phi0(y * hY, t * tau);
			U[t][Nx-1][y] = phi1(y * hY, t * tau);
			halftimeU[0][y] = phi0(y * hY, t * tau - tau / 2.0);
			halftimeU[Nx-1][y] = phi1(y * hY, t * tau - tau / 2.0);
		}
		// stage 1
		for (int y = 1; y < Ny - 1; ++y) {
			std::vector<double> lowD(Nx - 2);
			std::vector<double> mainD(Nx - 2);
			std::vector<double> upD(Nx - 2);
			std::vector<double> b(Nx - 2);

			mainD[0] = 2 * hX * hX * hY * hY + 2 * a * tau * hY * hY;
			upD[0] = -a * tau * hY * hY;
			for (int i = 1; i < Nx - 3; ++i) {
				lowD[i] = -a * tau * hY * hY;
				mainD[i] = 2 * hX * hX * hY * hY + 2 * a * tau * hY * hY;
				upD[i] = -a * tau * hY * hY;
			}
			lowD[Nx - 3] = -a * tau * hY * hY;
			mainD[Nx - 3] = 2 * hX * hX * hY * hY + 2 * a * tau * hY * hY;
			for (int x = 1; x < Nx - 1; ++x) {
				b[x - 1] =
					U[t - 1][x][y - 1] * a * tau * hX * hX
					+ U[t - 1][x][y] * (2 * hX * hX * hY * hY - 2 * a * tau * hX * hX)
					+ U[t - 1][x][y + 1] * a * tau * hX * hX;

			}
			b[0] -= (-a * tau * hY * hY) * phi0(y * hY, t * tau - tau / 2);
			b[Nx - 3] -= (-a * tau * hY * hY) * phi1(y * hY, t * tau - tau / 2);

			std::vector<double> curColumn(Nx - 2);
			tridiagonalAlgo(lowD, mainD, upD, b, curColumn, 0, 0, 0);
			for (int x = 1; x < Nx - 1; ++x) {
				halftimeU[x][y] = curColumn[x - 1];
			}
		}
		// stage 2
		for (int x = 1; x < Nx - 1; ++x) {
			std::vector<double> lowD(Ny - 2);
			std::vector<double> mainD(Ny - 2);
			std::vector<double> upD(Ny - 2);
			std::vector<double> b(Ny - 2);

			mainD[0] = 2 * hX * hX * hY * hY + 2 * a * tau * hX * hX;
			upD[0] = -a * tau * hX * hX;
			for (int i = 1; i < Ny - 3; ++i) {
				lowD[i] = -a * tau * hX * hX;
				mainD[i] = 2 * hX * hX * hY * hY + 2 * a * tau * hX * hX;
				upD[i] = -a * tau * hX * hX;
			}
			lowD[Ny - 3] = -a * tau * hX * hX;
			mainD[Ny - 3] = 2 * hX * hX * hY * hY + 2 * a * tau * hX * hX;
			for (int y = 1; y < Ny - 1; ++y) {
				b[y - 1] =
					halftimeU[x - 1][y] * a * tau * hY * hY
					+ halftimeU[x][y] * (2 * hX * hX * hY * hY - 2 * a * tau * hY * hY)
					+ halftimeU[x + 1][y] * a * tau * hY * hY;

			}
			b[0] -= (-a * tau * hX * hX) * phi2(x * hX, t * tau);
			b[Ny - 3] -= (-a * tau * hX * hX) * phi3(x * hX, t * tau);

			std::vector<double> curRow(Ny - 2);
			tridiagonalAlgo(lowD, mainD, upD, b, curRow, 0, 0, 0);
			for (int y = 1; y < Ny - 1; ++y) {
				U[t][x][y] = curRow[y - 1];
			}
		}

	}
	return U;
}
gridF fractSteps() {
	gridF U = std::vector<std::vector<std::vector<double>>>(Nt,
		std::vector<std::vector<double>>(Nx, std::vector<double>(Ny, 0.0)));

	for (int x = 0; x < Nx; ++x) {
		for (int y = 0; y < Ny; ++y) {
			U[0][x][y] = psi(x * hX, y * hY);
		}
	}

	for (int t = 1; t < Nt; ++t) {
		std::vector<std::vector<double>> halftimeU = std::vector<std::vector<double>>(Nx, std::vector<double>(Ny, 0.0));

		for (int x = 0; x < Nx; ++x) {
			U[t][x][0] = phi2(x * hX, t * tau);
			U[t][x][Ny - 1] = phi3(x * hX, t * tau);
			halftimeU[x][0] = phi2(x * hX, t * tau - tau / 2.0);
			halftimeU[x][Ny - 1] = phi3(x * hX, t * tau - tau / 2.0);
		}
		for (int y = 0; y < Ny; ++y) {
			U[t][0][y] = phi0(y * hY, t * tau);
			U[t][Nx - 1][y] = phi1(y * hY, t * tau);
			halftimeU[0][y] = phi0(y * hY, t * tau - tau / 2.0);
			halftimeU[Nx - 1][y] = phi1(y * hY, t * tau - tau / 2.0);
		}
		// stage 1
		for (int y = 1; y < Ny - 1; ++y) {
			std::vector<double> lowD(Nx - 2);
			std::vector<double> mainD(Nx - 2);
			std::vector<double> upD(Nx - 2);
			std::vector<double> b(Nx - 2);

			mainD[0] = hX * hX + 2 * a * tau;
			upD[0] = -a * tau;
			for (int i = 1; i < Nx - 3; ++i) {
				lowD[i] = -a * tau;
				mainD[i] = hX * hX + 2 * a * tau;
				upD[i] = -a * tau;
			}
			lowD[Nx - 3] = -a * tau;
			mainD[Nx - 3] = hX * hX + 2 * a * tau;
			for (int x = 1; x < Nx - 1; ++x) {
				b[x - 1] = U[t - 1][x][y] * hX * hX;

			}
			b[0] -= (-a * tau) * phi0(y * hY, t * tau - tau / 2);
			b[Nx - 3] -= (-a * tau) * phi1(y * hY, t * tau - tau / 2);

			std::vector<double> curColumn(Nx - 2);
			tridiagonalAlgo(lowD, mainD, upD, b, curColumn, 0, 0, 0);
			for (int x = 1; x < Nx - 1; ++x) {
				halftimeU[x][y] = curColumn[x - 1];
			}
		}
		// stage 2
		for (int x = 1; x < Nx - 1; ++x) {
			std::vector<double> lowD(Ny - 2);
			std::vector<double> mainD(Ny - 2);
			std::vector<double> upD(Ny - 2);
			std::vector<double> b(Ny - 2);

			mainD[0] = hY * hY + 2 * a * tau;
			upD[0] = -a * tau;
			for (int i = 1; i < Ny - 3; ++i) {
				lowD[i] = -a * tau;
				mainD[i] = hY * hY + 2 * a * tau;
				upD[i] = -a * tau;
			}
			lowD[Ny - 3] = -a * tau;
			mainD[Ny - 3] = hY * hY + 2 * a * tau;
			for (int y = 1; y < Ny - 1; ++y) {
				b[y - 1] = halftimeU[x][y] * hY * hY;

			}
			b[0] -= (-a * tau) * phi2(x * hX, t * tau);
			b[Ny - 3] -= (-a * tau) * phi3(x * hX, t * tau);

			std::vector<double> curRow(Ny - 2);
			tridiagonalAlgo(lowD, mainD, upD, b, curRow, 0, 0, 0);
			for (int y = 1; y < Ny - 1; ++y) {
				U[t][x][y] = curRow[y - 1];
			}
		}

	}
	return U;
}

// Error function
std::vector<double> getError(gridF& U) {
	std::vector<double> error(Nt, 0.0);
	for (int t = 0; t < Nt; ++t) {
		for (int x = 0; x < Nx; ++x) {
			for (int y = 0; y < Ny; ++y) {
				double curError = std::abs(U[t][x][y] - analSol(x * hX, y * hY, t * tau));
				if (curError > error[t]) {
					error[t] = curError;
				}
			}
		}
	}
	return error;
}
