#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

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
    ofstream outputFile("results_cn.txt");

    // Инициализация сетки и начальных условий
    double** u = new double*[Nx];
    for (int i = 0; i < Nx; ++i) {
        u[i] = new double[Nt];
    }

    // Инициализация начальных условий
    for (int i = 0; i < Nx; ++i) {
        double x = i * h;
        u[i][0] = sin(x);
    }

    // Схема Кранка-Николсона
    double alpha = a * k / (2 * h * h);
    double beta = c * k / 2;

    for (int j = 0; j < Nt - 1; ++j) {
        double t = (j + 1) * k;

        // Заполнение матрицы системы линейных уравнений
        double** A = new double*[Nx];
        for (int i = 0; i < Nx; ++i) {
            A[i] = new double[Nx];
        }

        for (int i = 0; i < Nx; ++i) {
            for (int l = 0; l < Nx; ++l) {
                A[i][l] = 0.0;
            }
            A[i][i] = 1 + 2 * alpha - beta;
            if (i > 0) {
                A[i][i - 1] = -alpha;
            }
            if (i < Nx - 1) {
                A[i][i + 1] = -alpha;
            }
        }

        // Заполнение вектора правой части
        double* b = new double[Nx];
        for (int i = 0; i < Nx; ++i) {
            b[i] = u[i][j] + alpha * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j])
                   + beta * (u[i + 1][j] + u[i - 1][j]);
        }

        // Решение системы линейных уравнений методом прогонки
        for (int i = 1; i < Nx; ++i) {
            double m = A[i][i - 1] / A[i - 1][i - 1];
            A[i][i] -= m * A[i - 1][i];
            b[i] -= m * b[i - 1];
        }

        u[Nx - 1][j + 1] = b[Nx - 1] / A[Nx - 1][Nx - 1];

        for (int i = Nx - 2; i >= 0; --i) {
            u[i][j + 1] = (b[i] - A[i][i + 1] * u[i + 1][j + 1]) / A[i][i];
        }

        // Очистка памяти
        for (int i = 0; i < Nx; ++i) {
            delete[] A[i];
        }
        delete[] A;
        delete[] b;
    }

    // Запись результатов в файл
    for (int j = 0; j < Nt; ++j) {
        double t = j * k;
        for (int i = 0; i < Nx; ++i) {
            double x = i * h;
            outputFile << x << " " << t << " " << u[i][j] << " " << analyticalSolution(x, t, a, c) << endl;
        }
    }

    // Очистка памяти
    for (int i = 0; i < Nx; ++i) {
        delete[] u[i];
    }
    delete[] u;

    // Закрытие файла
    outputFile.close();

    return 0;
}
