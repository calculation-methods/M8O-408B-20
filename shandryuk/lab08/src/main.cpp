#include <iostream>
#include <fstream>
#include <cmath>

const double Lx = 1.0; // Размеры области
const double Ly = 1.0;
const double T = 1.0;  // Время
const double alpha = 1.0;  // Коэффициент в уравнении

const int Nx = 50;  // Число узлов по x
const int Ny = 50;  // Число узлов по y
const int Nt = 100; // Число шагов по времени

// Функция для инициализации начальных условий
void initialize(double*** U, double*** V);

// Функция для решения начально-краевой задачи
void solve(double*** U, double*** V);

// Функция для вычисления аналитического решения
double analyticalSolution(double x, double y, double t);

// Функция для вычисления погрешности
double calculateError(double*** U, double*** V);

// Функция для записи результатов в файл
void writeToFile(double*** U, double*** V);

int main() {
    // Выделение памяти под массивы
    double*** U = new double**[Nx];
    double*** V = new double**[Nx];
    for (int i = 0; i < Nx; ++i) {
        U[i] = new double*[Ny];
        V[i] = new double*[Ny];
        for (int j = 0; j < Ny; ++j) {
            U[i][j] = new double[Nt + 1];
            V[i][j] = new double[Nt + 1];
        }
    }

    // Инициализация начальных условий
    initialize(U, V);

    // Решение начально-краевой задачи
    solve(U, V);

    // Запись результатов в файл
    writeToFile(U, V);

    // Освобождение выделенной памяти
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            delete[] U[i][j];
            delete[] V[i][j];
        }
        delete[] U[i];
        delete[] V[i];
    }
    delete[] U;
    delete[] V;

    return 0;
}

void initialize(double*** U, double*** V) {
    // Инициализация начальных условий
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            U[i][j][0] = i * Lx / (Nx - 1) * j * Ly / (Ny - 1); // Начальное условие
            V[i][j][0] = U[i][j][0];                             // Начальное условие
        }
    }
}

void solve(double*** U, double*** V) {
    // Решение начально-краевой задачи
    double hx = Lx / (Nx - 1);
    double hy = Ly / (Ny - 1);
    double tau = T / Nt;

    for (int k = 0; k < Nt; ++k) {
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                U[i][j][k + 1] = U[i][j][k] +
                    alpha * tau / (hx * hx) * (U[i + 1][j][k] - 2 * U[i][j][k] + U[i - 1][j][k]) +
                    alpha * tau / (hy * hy) * (U[i][j + 1][k] - 2 * U[i][j][k] + U[i][j - 1][k]) -
                    tau * i * hx * j * hy * sin((k + 1) * tau);
            }
        }
    }
}

double analyticalSolution(double x, double y, double t) {
    return x * y * cos(t);
}

double calculateError(double*** U, double*** V) {
    double error = 0.0;

    for (int k = 0; k <= Nt; ++k) {
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                double analytical = analyticalSolution(i * Lx / (Nx - 1), j * Ly / (Ny - 1), k * T / Nt);
                error += pow(U[i][j][k] - analytical, 2);
            }
        }
    }

    return sqrt(error / (Nx * Ny * (Nt + 1)));
}

void writeToFile(double*** U, double*** V) {
    // Запись результатов в файл
    std::ofstream outputFile("results.txt");
    for (int k = 0; k <= Nt; ++k) {
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                double analytical = analyticalSolution(i * Lx / (Nx - 1), j * Ly / (Ny - 1), k * T / Nt);
                outputFile << i * Lx / (Nx - 1) << " " << j * Ly / (Ny - 1) << " " << (k * T / Nt) << " " << U[i][j][k] << " " << analytical << " " << U[i][j][k] - analytical << std::endl;
            }
        }
    }
    outputFile.close();
}
