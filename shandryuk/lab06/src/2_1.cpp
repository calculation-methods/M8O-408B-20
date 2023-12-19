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
        U[0][j + 1] = 0;  // Двухточечная аппроксимация с первым порядком
        U[N][j + 1] = 0;  // Двухточечная аппроксимация с первым порядком
    }
}

// Решение системы линейных уравнений методом прогонки
void solveTridiagonalSystem(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, std::vector<double>& x) {
    int n = static_cast<int>(a.size());

    // Прямой ход
    for (int i = 1; i < n; ++i) {
        double m = a[i] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }

    // Обратный ход
    x[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
}

// Неявная схема с двухточечной аппроксимацией
void implicitScheme(std::vector<std::vector<double>>& U, double tau, double h, int N, int M) {
    double sigma = tau / h;

    // Решение системы на каждом временном шаге
    for (int j = 0; j < M; ++j) {
        // Создание трехдиагональных матриц
        std::vector<double> a(N - 1, -sigma / 2.0);
        std::vector<double> b(N - 1, 1 + sigma);
        std::vector<double> c(N - 1, -sigma / 2.0);
        std::vector<double> d(N - 1);

        // Заполнение вектора правой части
        for (int i = 1; i < N; ++i) {
            d[i - 1] = U[i][j] + tau * (2 * U[i][j] - 3 * U[i][j] +
                        exp(-i * h - (j + 1) * tau) * sin(i * h) * sin(2 * (j + 1) * tau));
        }

        // Решение системы линейных уравнений
        std::vector<double> x(N - 1);
        solveTridiagonalSystem(a, b, c, d, x);

        // Запись результатов в сетку
        for (int i = 1; i < N; ++i) {
            U[i][j + 1] = x[i - 1];
        }

        // Граничные условия
        U[0][j + 1] = 0;
        U[N][j + 1] = 0;
    }
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

    // Инициализация сетки
    std::vector<std::vector<double>> U(N + 1, std::vector<double>(M + 1, 0.0));

    // Инициализация начальных условий
    for (int i = 0; i <= N; ++i) {
        U[i][0] = 0.0;  // Начальное условие
    }

    // Решение с использованием явной схемы "крест"
    explicitCrossScheme(U, tau, h, N, M);

    // Вывод результатов и погрешности
    std::cout << "Results for Explicit Cross Scheme:" << std::endl;
    printResults(U, h, tau, N, M);

    // Очистка сетки
    U.assign(N + 1, std::vector<double>(M + 1, 0.0));

    // Решение с использованием неявной схемы
    implicitScheme(U, tau, h, N, M);

    // Вывод результатов и погрешности
    std::cout << "Results for Implicit Scheme:" << std::endl;
    printResults(U, h, tau, N, M);

    return 0;
}
