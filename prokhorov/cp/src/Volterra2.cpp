// Volterra2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

std::ostream& operator << (std::ostream& out, const std::vector<double>& vector) {
    for (double element : vector) {
        out << element << " ";
    }
    return out;
}

double a = 0, b = 1, n = 51, h = (b - a) / (n - 1);

double remainder(double x) {
    return x;
}

double kernel(double x, double s) {
    return sin(x - s) * 4. - x + s;
}

double solution(double x) {
    return (exp(x) + exp(-x)) * x / 2.;
}

std::vector<double> quadratureMethod(const std::vector<std::vector<double>>& K, const std::vector<double>& f, const std::vector<double>& A, 
    const std::vector<double>& x) {
    std::vector<double> answer(f);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            answer[i] += A[j] * K[i][j] * answer[j];
        }
        answer[i] /= (1 - A[i] * K[i][i]);
    }
    return answer;
}

double norm(const std::vector<double>& x) {
    double sum = 0.;
    for (double number: x) {
        sum += number * number;
    }
    return sqrt(sum);
}

const std::vector<double> operator - (const std::vector<double>& x, const std::vector<double>& y) {
    std::vector<double> answer(n);
    for (int i = 0; i < n; ++i) {
        answer[i] = x[i] - y[i];
    }
    return answer;
}

std::vector<double> iterationMethod(const std::vector<std::vector<double>>& K, const std::vector<double>& f, const std::vector<double>& A,
    const std::vector<double>& x, double epsilon = 1e-9) {
    std::vector<double> next(f);
    while (true) {
        std::vector<double> current(next);
        for (int i = 0; i < n; ++i) {
            next[i] = f[i];
            for (int j = 0; j < i + 1; ++j) {
                next[i] += A[j] * K[i][j] * current[j];
            }
        }
        if (norm(next - current) / norm(next) < epsilon) {
            break;
        }
    }
    return next;
}

void output_to_file(std::string path, const std::vector<double>& vector) {
    std::ofstream fout(path);
    for (double element : vector) {
        fout << element << " ";
    }
    fout << "\n";
    fout.close();
}

int main()
{
    std::setlocale(LC_ALL, "Russian");
    std::vector<double> trapezoid(n, h), Simpson(n, h / 3.), x(n, a), f(n);
    trapezoid[0] /= 2;
    trapezoid[n - 1] /= 2;
    f[0] = remainder(x[0]);
    for (int i = 1; i < n; ++i) {
        x[i] = x[i - 1] + h;
        f[i] = remainder(x[i]);
        if (i < n - 1) {
            Simpson[i] *= 2;
            if (i % 2) {
                Simpson[i] *= 2;
            }
        }
    }
    std::vector<std::vector<double>> K(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i + 1; ++j) {
            K[i][j] = kernel(x[i], x[j]);
        }
    }
    std::vector<double> answer(n);
    for (int i = 0; i < n; ++i) {
        answer[i] = solution(x[i]);
    }
    std::cout << "1) Истинное решение: " << answer << "\n";
    output_to_file("analytical_solution.txt", answer);
    std::vector<double> quadratureTrapezoidSolution = quadratureMethod(K, f, trapezoid, x);
    std::cout << "2) Решение методом квадратур с использованием метода трапеций: " << quadratureTrapezoidSolution << "\n";
    std::cout << "Квадратичная норма между решениями: " << norm(quadratureTrapezoidSolution - answer) << "\n";
    output_to_file("quadrature_trapezoid.txt", quadratureTrapezoidSolution);
    std::vector<double> quadratureSimpsonSolution = quadratureMethod(K, f, Simpson, x);
    std::cout << "3) Решение методом квадратур с использованием метода Симпсона: " << quadratureSimpsonSolution << "\n";
    std::cout << "Квадратичная норма между решениями: " << norm(quadratureSimpsonSolution - answer) << "\n";
    output_to_file("quadrature_Simpson.txt", quadratureSimpsonSolution);
    std::vector<double> iterationTrapezoidSolution = iterationMethod(K, f, trapezoid, x);
    std::cout << "4) Решение методом простых итераций с использованием метода трапеций: " << iterationTrapezoidSolution << "\n";
    std::cout << "Квадратичная норма между решениями: " << norm(iterationTrapezoidSolution - answer) << "\n";
    output_to_file("iterations_trapezoid.txt", iterationTrapezoidSolution);
    std::vector<double> iterationSimpsonSolution = iterationMethod(K, f, Simpson, x);
    std::cout << "5) Решение методом простых итераций с использованием метода Симпсона: " << iterationSimpsonSolution << "\n";
    std::cout << "Квадратичная норма между решениями: " << norm(iterationSimpsonSolution - answer) << "\n";
    output_to_file("iterations_Simpson.txt", iterationSimpsonSolution);
    return 0;
}
