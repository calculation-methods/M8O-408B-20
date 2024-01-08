#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>

class SparseMatrix {
public:
    SparseMatrix(int rows, int cols) : rows(rows), cols(cols) {}

    void addValue(int row, int col, double value) {
        elements.emplace_back(row, col, value);
    }

    void setB(const std::vector<double>& b_values) {
        b = b_values;
    }

    int getRows() const {
        return rows;
    }

    int getCols() const {
        return cols;
    }

    const std::vector<std::tuple<int, int, double>>& getElements() const {
        return elements;
    }

    const std::vector<double>& getB() const {
        return b;
    }

private:
    int rows;
    int cols;
    std::vector<std::tuple<int, int, double>> elements;
    std::vector<double> b;
};

class Solver {
public:
    Solver(const SparseMatrix& matrix, const std::vector<double>& b,
           const std::string& output_file, const std::vector<double>& x0 = {},
           double eps = 1e-5)
        : output(output_file.empty() ? "res_default" : output_file),
          matrix(matrix), b(b), eps(eps), shape(matrix.getRows()), x0(x0), k(0) {}

    std::vector<double> solve(int max_iter = 100000) {
        std::vector<double> x_new, r_new, p_new;
        std::vector<double> x0 = this->x0;
        std::vector<double> r0 = vectorSubtract(b, matrixMultiply(matrix, x0));
        std::vector<double> p0 = r0;

        for (int iter = 0; iter < max_iter; ++iter) {
            std::vector<double> temp = matrixMultiply(matrix, p0);
            double norm_0 = vectorDot(r0, r0);
            double alpha_i = norm_0 / vectorDot(temp, p0);
            x_new = vectorAdd(x0, vectorMultiply(p0, alpha_i));
            r_new = vectorSubtract(r0, vectorMultiply(temp, alpha_i));
            double norm_new = vectorDot(r_new, r_new);
            double beta_i = norm_new / norm_0;
            p_new = vectorAdd(r_new, vectorMultiply(p0, beta_i));

            r0 = r_new;
            p0 = p_new;
            x0 = x_new;

            k += 1;

            if (vectorNorm(r_new) < eps) {
                break;
            }
        }
        return x0;
    }

    void solveAndPrint() {
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<double> x = solve();
        auto end = std::chrono::high_resolution_clock::now();
        auto start2 = std::chrono::high_resolution_clock::now();
        std::vector<double> x2 = solveLinearSystem(matrix, b);
        auto end2 = std::chrono::high_resolution_clock::now();

        std::cout << "Custom solution:\n";
        printVector(x);
        std::cout << "eps=" << eps << " shape=" << shape << " iterations=" << k
                  << " mean=" << vectorMean(x) << " time=" << std::chrono::duration<double>(end - start).count() << " seconds\n";

        std::cout << "NumPy solution:\n";
        printVector(x2);
        std::cout << "mean=" << vectorMean(x2) << " time=" << std::chrono::duration<double>(end2 - start2).count() << " seconds\n";
    }

private:
    std::string output;
    SparseMatrix matrix;
    std::vector<double> b;
    double eps;
    int shape;
    std::vector<double> x0;
    int k;

    static std::vector<double> vectorAdd(const std::vector<double>& a, const std::vector<double>& b) {
        std::vector<double> result(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    static std::vector<double> vectorSubtract(const std::vector<double>& a, const std::vector<double>& b) {
        std::vector<double> result(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    static std::vector<double> vectorMultiply(const std::vector<double>& a, double scalar) {
        std::vector<double> result(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            result[i] = a[i] * scalar;
        }
        return result;
    }

    static double vectorDot(const std::vector<double>& a, const std::vector<double>& b) {
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            result += a[i] * b[i];
        }
        return result;
    }

    static double vectorNorm(const std::vector<double>& v) {
        double result = 0.0;
        for (double val : v) {
            result += val * val;
        }
        return std::sqrt(result);
    }

    static double vectorMean(const std::vector<double>& v) {
        double sum = 0.0;
        for (double val : v) {
            sum += val;
        }
        return sum / v.size();
    }

    static void printVector(const std::vector<double>& v) {
        for (double val : v) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }

    static std::vector<double> solveLinearSystem(const SparseMatrix& matrix, const std::vector<double>& b) {
        // Simple Gaussian elimination method for solving Ax = b.
        int n = matrix.getRows();
        std::vector<std::vector<double>> augmentedMatrix(n, std::vector<double>(n + 1));

        for (const auto& element : matrix.getElements()) {
            int row = std::get<0>(element);
            int col = std::get<1>(element);
            double value = std::get<2>(element);
            augmentedMatrix[row][col] = value;
        }

        for (int i = 0; i < n; ++i) {
            augmentedMatrix[i][n] = b[i];
        }

        for (int i = 0; i < n; ++i) {
            // Pivot for the current column
            double pivot = augmentedMatrix[i][i];

            // Make the diagonal element 1
            for (int j = i + 1; j <= n; ++j) {
                augmentedMatrix[i][j] /= pivot;
            }

            for (int k = 0; k < n; ++k) {
                if (k != i) {
                    double factor = augmentedMatrix[k][i];
                    for (int j = i; j <= n; ++j) {
                        augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                    }
                }
            }
        }

        std::vector<double> solution(n);
        for (int i = 0; i < n; ++i) {
            solution[i] = augmentedMatrix[i][n];
        }

        return solution;
    }

    static std::vector<double> matrixMultiply(const SparseMatrix& matrix, const std::vector<double>& v) {
        std::vector<double> result(matrix.getRows(), 0.0);
        for (const auto& element : matrix.getElements()) {
            int row = std::get<0>(element);
            int col = std::get<1>(element);
            double value = std::get<2>(element);
            result[row] += value * v[col];
        }
        return result;
    }
};

SparseMatrix readMatrix(const std::string& filename) {
    std::ifstream file(filename);
    int shape;
    file >> shape;

    SparseMatrix matrix(shape, shape);

    for (int i = 0; i < shape; ++i) {
        for (int j = 0; j < shape; ++j) {
            double value;
            file >> value;
            if (value != 0.0) {
                matrix.addValue(i, j, value);
            }
        }
    }

    std::vector<double> b(shape);
    for (int i = 0; i < shape; ++i) {
        file >> b[i];
    }

    matrix.setB(b);

    return matrix;
}

int main() {
    std::string input_file = "matrix.txt";
    std::string output_file = "output_file.txt";
    double eps = 0.01;

    SparseMatrix matrix = readMatrix(input_file);
    std::vector<double> b = matrix.getB();

    Solver solver(matrix, b, output_file, {}, eps);
    solver.solveAndPrint();

    return 0;
}
