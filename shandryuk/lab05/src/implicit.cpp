#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const double pi = 3.141592653589793;

// Функция для вычисления аналитического решения
double analyticalSolution(double x, double t, double a, double c) {
    return exp((c - a) * t) * sin(x);
}

int main() {
    // Параметры задачи
    double a = 1.0;
    double c = 1.0;
    double L = pi / 2.0;
    double T = 1.0;

    // Шаги сетки
    double h = 0.1;  // шаг по x
    double k = 0.01; // шаг по t

    // Количество узлов сетки
    int Nx = static_cast<int>(L / h) + 1;
    int Nt = static_cast<int>(T / k) + 1;

    // Открытие файла для записи результатов
    ofstream outputFile("results_implicit.txt");

    // Инициализация сетки и начальных условий
    MatrixXd u(Nx, Nt);

    // Инициализация начальных условий
    for (int i = 0; i < Nx; ++i) {
        double x = i * h;
        u(i, 0) = sin(x);
    }

    // Неявная конечно-разностная схема
    double alpha = a * k / (h * h);
    double beta = c * k;

    // Создание трехдиагональной матрицы для решения системы линейных уравнений
    MatrixXd A(Nx, Nx);
    A.setZero();
    A.diagonal().setConstant(1 + 2 * alpha - beta);
    A.diagonal(1).setConstant(-alpha);
    A.diagonal(-1).setConstant(-alpha);

    for (int j = 0; j < Nt - 1; ++j) {
        VectorXd b(Nx);
        b.setZero();

        double t = (j + 1) * k;

        // Заполнение вектора правой части
        for (int i = 1; i < Nx - 1; ++i) {
            b(i) = u(i, j) + beta * h * h * analyticalSolution(i * h, t, a, c);
        }

        // Решение системы линейных уравнений Ax = b
        VectorXd x = A.fullPivLu().solve(b);

        // Запись результатов
        for (int i = 1; i < Nx - 1; ++i) {
            u(i, j + 1) = x(i);
        }
    }

    // Запись результатов в файл
    for (int j = 0; j < Nt; ++j) {
        double t = j * k;
        for (int i = 0; i < Nx; ++i) {
            double x = i * h;
            outputFile << x << " " << t << " " << u(i, j) << " " << analyticalSolution(x, t, a, c) << endl;
        }
    }

    // Закрытие файла
    outputFile.close();

    return 0;
}
