#include <iostream>
#include <cmath>
#include <vector>

// Аналитическое решение
double analyticalSolution(double x, double t) {
    return exp(-x - t) * sin(x) * sin(2 * t);
}

// Явная схема "крест" с двухточечной аппроксимацией
void explicitCrossScheme(std::vector<std::vector<double>>& U, double tau, double h, int N, int M) {
    double sigma = tau / h;

    for (int j = 0; j < M; ++j) {
        for (int i = 1; i < N; ++i) {
            U[i][j + 1] = sigma * (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) +
                           2 * U[i][j] - U[i][j - 1] + tau * (2 * U[i][j] - 3 * U[i][j] +
                           exp(-i * h - (j + 1) * tau) * sin(i * h) * sin(2 * (j + 1) * tau));
        }

        // Граничные условия
        U[0][j + 1] = 0;
        U[N][j + 1] = 0;
    }
}

// Неявная схема с двухточечной аппроксимацией
void implicitScheme(std::vector<std::vector<double>>& U, double tau, double h, int N, int M) {
    double sigma = tau / h;

    for (int j = 0; j < M; ++j) {
        for (int i = 1; i < N; ++i) {
            U[i][j + 1] = (2 * U[i][j] - U[i][j - 1] + tau * (2 * U[i][j] - 3 * U[i][j] +
                           exp(-i * h - (j + 1) * tau) * sin(i * h) * sin(2 * (j + 1) * tau))) /
                           (1 + 2 * sigma);
        }

        // Граничные условия
        U[0][j + 1] = 0;
        U[N][j + 1] = 0;
    }
}

// Функция для вычисления погрешности
double calculateError(const std::vector<std::vector<double>>& numerical, const std::vector<std::vector<double>>& analytical) {
    double maxError = 0.0;
    for (size_t i = 0; i < numerical.size(); ++i) {
        for (size_t j = 0; j < numerical[i].size(); ++j) {
            double error = std::abs(numerical[i][j] - analyticalSolution(i, j));
            maxError = std::max(maxError, error);
        }
    }
    return maxError;
}

// Функция для вывода результатов
void printResults(const std::vector<std::vector<double>>& U, double h, double tau, int N, int M) {
    for (int j = 0; j <= M; ++j) {
        double t = j * tau;
        for (int i = 0; i <= N; ++i) {
            double x = i * h;
            double analytical = analyticalSolution(x, t);
            double error = std::abs(U[i][j] - analytical);
            std::cout << "Time: " << t << ", Position: " << x << ", Numerical: " << U[i][j]
                      << ", Analytical: " << analytical << ", Error: " << error << std::endl;
        }
    }
}

int main() {
    // Параметры сетки
    double L = M_PI;
    int N = 100;
    int M = 100;
    double T = 1.0;
    double tau, h;

    // Выбор шагов
    std::cout << "Enter tau (time step): ";
    std::cin >> tau;
    std::cout << "Enter h (space step): ";
    std::cin >> h;

    // Количество временных и пространственных узлов
    M = static_cast<int>(T / tau);
    N = static_cast<int>(L / h);

    // Инициализация сеток
    std::vector<std::vector<double>> U_explicit(N + 1, std::vector<double>(M + 1, 0.0));
    std::vector<std::vector<double>> U_implicit(N + 1, std::vector<double>(M + 1, 0.0));

    // Инициализация начальных условий
    for (int i = 0; i <= N; ++i) {
        U_explicit[i][0] = 0.0;  // Начальное условие
        U_implicit[i][0] = 0.0;  // Начальное условие
    }

    // Решение с использованием явной схемы
    explicitCrossScheme(U_explicit, tau, h, N, M);

    // Решение с использованием неявной схемы
    implicitScheme(U_implicit, tau, h, N, M);

    // Вывод результатов и погрешности
    std::cout << "Results for Explicit Scheme:" << std::endl;
    printResults(U_explicit, h, tau, N, M);
    std::cout << "Max Error: " << calculateError(U_explicit, U_implicit) << std::endl;

    std::cout << "Results for Implicit Scheme:" << std::endl;
    printResults(U_implicit, h, tau, N, M);

    return 0;
}
