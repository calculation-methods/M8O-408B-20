// 7.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <cmath>

double phi0(double y) {
    return cos(y);
}

double phi1(double y) {
    return exp(1.) * cos(y);
}

double psi0(double x) {
    return exp(x);
}

double psi1(double x) {
    return 0.;
}

double solution(double x, double y) {
    return exp(x) * cos(y);
}

double L2(std::vector<double>& x, std::vector<double>& y) {
    int n = x.size();
    double l2 = 0.;
    for (int i = 0; i < n; ++i) {
        l2 += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sqrt(l2);
}

std::vector<double> multiplication(std::vector<std::vector<double>>& A, std::vector<double>& b) {
    int n = A.size();
    std::vector<double> answer(n, 0.);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            answer[i] += A[i][j] * b[j];
        }
    }
    return answer;
}

void addition(std::vector<double>& x, std::vector<double>& y) {
    int n = x.size();
    for (int i = 0; i < n; ++i) {
        x[i] += y[i];
    }
}

std::vector<std::vector<double>> analytical_solution(double x_begin, double x_end, double y_begin,
    double y_end, double hx, double hy) {
    int length_x = static_cast<int>((x_end - x_begin) / hx) + 1;
    int length_y = static_cast<int>((y_end - y_begin) / hy) + 1;
    std::vector<double> x(length_x), y(length_y);
    x[0] = x_begin;
    for (int i = 1; i < length_x; ++i) {
        x[i] = x[i - 1] + hx;
    }
    y[0] = y_begin;
    for (int i = 1; i < length_y; ++i) {
        y[i] = y[i - 1] + hy;
    }
    std::vector<std::vector<double>> result(length_x, std::vector<double>(length_y));
    for (int i = 0; i < length_x; ++i) {
        for (int j = 0; j < length_y; ++j) {
            result[i][j] = solution(x[i], y[j]);
        }
    }
    return result;
}

std::pair<std::vector<std::vector<double>>, int> finite_difference_method(double x_begin, double x_end, double y_begin,
    double y_end, double hx, double hy, double epsilon,
    std::function<std::pair<std::vector<double>, int>(std::vector<std::vector<double>>&, std::vector<double>&, double)> method) {
    int length_x = static_cast<int>((x_end - x_begin) / hx) + 1;
    int length_y = static_cast<int>((y_end - y_begin) / hy) + 1;
    std::vector<double> x(length_x), y(length_y);
    x[0] = x_begin;
    for (int i = 1; i < length_x; ++i) {
        x[i] = x[i - 1] + hx;
    }
    y[0] = y_begin;
    for (int i = 1; i < length_y; ++i) {
        y[i] = y[i - 1] + hy;
    }
    std::vector<std::vector<double>> result(length_x, std::vector<double>(length_y));
    for (int i = 0; i < length_x; ++i) {
        result[i][0] = psi0(x[i]);
        result[i][length_y - 1] = psi1(x[i]);
    }
    for (int i = 0; i < length_y; ++i) {
        result[0][i] = phi0(y[i]);
        result[length_x - 1][i] = phi1(y[i]);
    }
    std::vector<std::vector<int>> mapping(length_x, std::vector<int>(length_y));
    int current_equation = 0;
    for (int i = 1; i < length_x - 1; ++i) {
        for (int j = 1; j < length_y - 1; ++j) {
            mapping[i][j] = current_equation++;
        }
    }
    int number_of_equations = (length_x - 2) * (length_y - 2);
    std::vector<std::vector<double>> A(number_of_equations, std::vector<double>(number_of_equations, 0.));
    std::vector<double> b(number_of_equations, 0.);
    for (int i = 1; i < length_x - 1; ++i) {
        for (int j = 1; j < length_y - 1; ++j) {
            current_equation = mapping[i][j];
            A[current_equation][mapping[i][j]] = 1;
            if (j == 1) {
                b[current_equation] += psi0(x[i]) * hx * hx / ((hx * hx + hy * hy) * 2.);
            }
            else {
                A[current_equation][mapping[i][j - 1]] = -hx * hx / ((hx * hx + hy * hy) * 2.);
            }
            if (j == length_y - 2) {
                b[current_equation] += psi1(x[i]) * hx * hx / ((hx * hx + hy * hy) * 2.);
            }
            else {
                A[current_equation][mapping[i][j + 1]] = -hx * hx / ((hx * hx + hy * hy) * 2.);
            }
            if (i == 1) {
                b[current_equation] += phi0(y[j]) * hy * hy / ((hx * hx + hy * hy) * 2.);
            }
            else {
                A[current_equation][mapping[i - 1][j]] = -hy * hy / ((hx * hx + hy * hy) * 2.);
            }
            if (i == length_x - 2) {
                b[current_equation] += phi1(y[j]) * hy * hy / ((hx * hx + hy * hy) * 2.);
            }
            else {
                A[current_equation][mapping[i + 1][j]] = -hy * hy / ((hx * hx + hy * hy) * 2.);
            }
        }
    }
    auto [answer, iterations] = method(A, b, epsilon);
    for (int i = 1; i < length_x - 1; ++i) {
        for (int j = 1; j < length_y - 1; ++j) {
            result[i][j] = answer[mapping[i][j]];
        }
    }
    return std::make_pair(result, iterations);
}

std::pair<std::vector<double>, int> iterative(std::vector<std::vector<double>>& A, std::vector<double>& b, double epsilon) {
    int n = A.size();
    std::vector<std::vector<double>> alpha(n, std::vector<double>(n, 0.));
    std::vector<double> beta(n, 0.);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                alpha[i][j] -= A[i][j] / A[i][i];
            }
        }
        beta[i] = b[i] / A[i][i];
    }
    int iterations = 0;
    std::vector<double> current_x(beta);
    bool converge = false;
    while (!converge) {
        std::vector<double> previous_x(current_x);
        current_x = multiplication(alpha, previous_x);
        addition(current_x, beta);
        ++iterations;
        double l2 = L2(previous_x, current_x);
        converge = l2 <= epsilon;
        std::cout << "Iterative method. Iteration: " << iterations <<  ", l2: " << l2 << "\n";
    }
    return std::make_pair(current_x, iterations);
}

std::vector<double> seidel_multiplication(std::vector<std::vector<double>>& alpha, std::vector<double>& x, std::vector<double>& beta) {
    int n = alpha.size(), m = alpha[0].size();
    std::vector<double> result(x);
    for (int i = 0; i < n; ++i) {
        result[i] = beta[i];
        for (int j = 0; j < m; ++j) {
            result[i] += alpha[i][j] * result[j];
        }
    }
    return result;
}

std::pair<std::vector<double>, int> seidel(std::vector<std::vector<double>>& A, std::vector<double>& b, double epsilon) {
    int n = A.size();
    std::vector<std::vector<double>> alpha(n, std::vector<double>(n, 0.));
    std::vector<double> beta(n, 0.);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                alpha[i][j] -= A[i][j] / A[i][i];
            }
        }
        beta[i] = b[i] / A[i][i];
    }
    int iterations = 0;
    std::vector<double> current_x(beta);
    bool converge = false;
    while (!converge) {
        std::vector<double> previous_x(current_x);
        current_x = seidel_multiplication(alpha, previous_x, beta);
        ++iterations;
        double l2 = L2(previous_x, current_x);
        converge = l2 <= epsilon;
        std::cout << "Seidel method. Iteration: " << iterations << ", l2: " << l2 << "\n";
    }
    return std::make_pair(current_x, iterations);
}

double w = 1.5;

void relaxations_additon(std::vector<double>& x, std::vector<double>& y) {
    int n = x.size();
    for (int i = 0; i < n; ++i) {
        x[i] = x[i] * w + y[i] * (1. - w);
    }
}

std::pair<std::vector<double>, int> relaxations(std::vector<std::vector<double>>& A, std::vector<double>& b, double epsilon) {
    int n = A.size();
    std::vector<std::vector<double>> alpha(n, std::vector<double>(n, 0.));
    std::vector<double> beta(n, 0.);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                alpha[i][j] -= A[i][j] / A[i][i];
            }
        }
        beta[i] = b[i] / A[i][i];
    }
    int iterations = 0;
    std::vector<double> current_x(beta);
    bool converge = false;
    while (!converge) {
        std::vector<double> previous_x(current_x);
        current_x = seidel_multiplication(alpha, previous_x, beta);
        relaxations_additon(current_x, previous_x);
        ++iterations;
        double l2 = L2(previous_x, current_x);
        converge = l2 <= epsilon;
        std::cout << "Relaxations method. Iteration: " << iterations << ", l2: " << l2 << "\n";
    }
    return std::make_pair(current_x, iterations);
}

double max_abs_error(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B) {
    int n = A.size(), m = A[0].size();
    double max = 0.;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            max = std::max(max, std::abs(A[i][j] - B[i][j]));
        }
    }
    return max;
}

double mean_abs_error(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B) {
    int n = A.size(), m = A[0].size();
    double mean = 0., prod = static_cast<double>(n * m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            mean += std::abs(A[i][j] - B[i][j]) / prod;
        }
    }
    return mean;
}

void output_to_file(std::string filepath, const std::vector<std::vector<double>>& arr) {
    std::ofstream fout(filepath);
    int n = arr.size(), m = arr[0].size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            fout << arr[i][j] << " ";
        }
        fout << "\n";
    }
    fout.close();
}

int main()
{
    double x_begin = 0., x_end = 1., y_begin = 0., y_end = acos(-1) / 2., hx = 0.01, hy = 0.01, epsilon = 1e-3;
    std::vector<std::vector<double>> as = analytical_solution(x_begin, x_end, y_begin, y_end, hx, hy);
    output_to_file("analytical_solution.txt", as);
    auto [is, ii] = finite_difference_method(x_begin, x_end, y_begin, y_end, hx, hy, epsilon, iterative);
    std::cout << "Iterative method took " << ii << " iterations until convergence\n";
    std::cout << "Max abs error between analytical solution and iterative method solution: " << max_abs_error(as, is) << "\n";
    std::cout << "Mean abs error between analytical solution and iterative method solution: " << mean_abs_error(as, is) << "\n";
    output_to_file("iterative_method.txt", is);
    auto [ss, si] = finite_difference_method(x_begin, x_end, y_begin, y_end, hx, hy, epsilon, seidel);
    std::cout << "Seidel method took " << si << " iterations until convergence\n";
    std::cout << "Max abs error between analytical solution and Seidel method solution: " << max_abs_error(as, ss) << "\n";
    std::cout << "Mean abs error between analytical solution and Seidel method solution: " << mean_abs_error(as, ss) << "\n";
    output_to_file( "seidel_method.txt", ss);
    auto [rs, ri] = finite_difference_method(x_begin, x_end, y_begin, y_end, hx, hy, epsilon, relaxations);
    std::cout << "Relaxations method took " << ri << " iterations until convergence\n";
    std::cout << "Max abs error between analytical solution and relaxations method solution: " << max_abs_error(as, rs) << "\n";
    std::cout << "Mean abs error between analytical solution and relaxations method solution: " << mean_abs_error(as, rs) << "\n";
    output_to_file("relaxations_method.txt", rs);
    return 0;
}


