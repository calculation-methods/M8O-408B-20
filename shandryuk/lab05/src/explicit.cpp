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
    ofstream outputFile("results_explicit.txt");

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

    // Явная конечно-разностная схема
    for (int j = 0; j < Nt - 1; ++j) {
        double t = j * k;
        for (int i = 1; i < Nx - 1; ++i) {
            double x = i * h;

            // Шаг схемы
            u[i][j + 1] = u[i][j] + a * k / (h * h) * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j])
                          + c * k * u[i][j];
        }
    }

    // Запись результатов в файл
    for (int j = 0; j < Nt; ++j) {
        double t = j * k;
        for (int i = 0; i < Nx; ++i) {
            double x = i * h;
            outputFile << x << " " << t << " " << u[i][j] << " " << analyticalSolution(x, t, a, c) << endl;
        }
    }

    // Закрытие файла
    outputFile.close();

    // Очистка памяти
    for (int i = 0; i < Nx; ++i) {
        delete[] u[i];
    }
    delete[] u;

    return 0;
}
