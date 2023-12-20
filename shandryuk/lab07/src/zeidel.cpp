#include <iostream>
#include <fstream>
#include <cmath>

// Глобальные параметры
const int Nx = 100; // Количество шагов по x
const int Ny = 100; // Количество шагов по y
const double Lx = M_PI / 2.0; // Граница по x
const double Ly = M_PI / 2.0; // Граница по y
const double hx = Lx / (Nx - 1); // Шаг по x
const double hy = Ly / (Ny - 1); // Шаг по y
const double epsilon = 1e-6; // Критерий останова метода Зейделя

// Аналитическое решение
double analyticalSolution(double x, double y) {
    return exp(-x) * cos(x) * cos(y);
}

int main() {
    // Инициализация сетки
    double u[Nx][Ny];
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            u[i][j] = 0.0;
        }
    }

    // Задание граничных условий
    for (int j = 0; j < Ny; ++j) {
        u[0][j] = cos(j * hy);
        u[Nx - 1][j] = 0.0;
    }
    for (int i = 0; i < Nx; ++i) {
        u[i][0] = exp(-i * hx) * cos(i * hx);
        u[i][Ny - 1] = 0.0;
    }

    // Метод Зейделя
    double error = epsilon + 1.0;
    int iteration = 0;
    while (error > epsilon) {
        error = 0.0;
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double newU = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]
                                      + hx * hy * (-2 * u[i][j] + 2 * cos(i * hx) + 3 * u[i][j]));
                error = std::max(error, std::abs(newU - u[i][j]));
                u[i][j] = newU;
            }
        }
        ++iteration;
    }

    // Вычисление погрешности и вывод результатов в файл
    std::ofstream outputFile("result_zeidel.txt");
    outputFile << "x\ty\tNumerical\tAnalytical\tError\n";
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double x = i * hx;
            double y = j * hy;
            double numericalSolution = u[i][j];
            double exactSolution = analyticalSolution(x, y);
            double currentError = std::abs(numericalSolution - exactSolution);
            outputFile << x << "\t" << y << "\t" << numericalSolution << "\t" << exactSolution << "\t" << currentError << "\n";
        }
    }
    outputFile.close();

    std::cout << "Calculation completed in " << iteration << " iterations." << std::endl;

    return 0;
}
