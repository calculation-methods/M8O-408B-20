#include <iostream>
#include <vector>
#include <cmath>
#define EPS 1e-7


// boundary conditions by x
double phi_0(double y) {
    return y;
}

double phi_1(double y) {
    return 1 + y;
}

// boundary conditions by y
double psi_0(double x) {
    return x;
}

double psi_1(double x) {
    return 1 + x;
}

// analytical solution
double solution(double x, double y) {
    return x + y;
}

// analytical solution
double solution(double x, double y) {
    return x + y;
}


std::vector<std::vector<double>> get_analytical_solution(std::pair<double, double> x_range, std::pair<double, double> y_range, double h_x, double h_y) {
    std::vector<double> x;
    std::vector<double> y;
    
    for (double i = x_range.first; i < x_range.second; i += h_x) {
        x.push_back(i);
    }
    
    for (double j = y_range.first; j < y_range.second; j += h_y) {
        y.push_back(j);
    }
    
    std::vector<std::vector<double>> res(x.size(), std::vector<double>(y.size()));
    
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < y.size(); j++) {
            res[i][j] = solution(x[i], y[j]);
        }
    }
    
    return res;
}

std::map<std::string, std::vector<std::vector<double>>> solutions;

double max_abs_error(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {    
    assert(A.size() == B.size());
    assert(A[0].size() == B[0].size());
    
    double max_error = 0.0;
    
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            double error = std::abs(A[i][j] - B[i][j]);
            if (error > max_error) {
                max_error = error;
            }
        }
    }
    
    return max_error;
}

double mean_abs_error(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    assert(A.size() == B.size());
    assert(A[0].size() == B[0].size());
    
    double sum_error = 0.0;
    
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            sum_error += std::abs(A[i][j] - B[i][j]);
        }
    }
    
    return sum_error / (A.size() * A[0].size());
}

double L2_norm(const std::vector<double>& X) {
    double l2_norm = 0.0;
    
    for (int i = 0; i < X.size(); i++) {
        l2_norm += X[i] * X[i];
    }
    
    return std::sqrt(l2_norm);
}


std::vector<std::vector<double>> finite_difference_schema(std::pair<double, double> x_range, std::pair<double, double> y_range, double h_x, double h_y, string method, double phi_0, double phi_1, double psi_0, double psi_1, double eps) {
    std::vector<double> x;
    std::vector<double> y;
    
    for (double i = x_range.first; i < x_range.second; i += h_x) {
        x.push_back(i);
    }
    
    for (double j = y_range.first; j < y_range.second; j += h_y) {
        y.push_back(j);
    }
    
    std::vector<std::vector<double>> res(x.size(), std::vector<double>(y.size()));
    
    // Step 1. Initialise grid with border conditions
    for (int cur_x_id = 0; cur_x_id < x.size(); cur_x_id++) {
        res[cur_x_id][0] = psi_0;
        res[cur_x_id][y.size() - 1] = psi_1;
    }
    
    for (int cur_y_id = 0; cur_y_id < y.size(); cur_y_id++) {
        res[0][cur_y_id] = phi_0(y[cur_y_id]);
        res[x.size() - 1][cur_y_id] = phi_1(y[cur_y_id]);
    }
    
    // Step 2. Create system of equations
    std::vector<std::vector<int>> mapping(x.size(), std::vector<int>(y.size(), 0));
    int cur_eq_id = 0;
    
    for (int cur_x_id = 1; cur_x_id < x.size() - 1; cur_x_id++) {
        for (int cur_y_id = 1; cur_y_id < y.size() - 1; cur_y_id++) {
            mapping[cur_x_id][cur_y_id] = cur_eq_id;
            cur_eq_id++;
        }
    }
    
    int nums_of_equations = (x.size() - 2) * (y.size() - 2);
    std::vector<std::vector<double>> A(nums_of_equations, std::vector<double>(nums_of_equations, 0.0));
    std::vector<double> b(nums_of_equations, 0.0);
    
    for (int cur_x_id = 1; cur_x_id < x.size() - 1; cur_x_id++) {
        for (int cur_y_id = 1; cur_y_id < y.size() - 1; cur_y_id++) {
            int cur_eq_id = mapping[cur_x_id][cur_y_id];
            A[cur_eq_id][cur_eq_id] = 1.0; // u_{i, j}
            
            if (cur_y_id - 1 == 0) {
                // u_{i, j-1} is already known from border conditions -> move the result to b
                b[cur_eq_id] += psi_0 * std::pow(h_x, 2) / (2 * (std::pow(h_x, 2) + std::pow(h_y, 2)));
            } else {
                A[cur_eq_id][mapping[cur_x_id][cur_y_id - 1]] = -std::pow(h_x, 2) / (2 * (std::pow(h_x, 2) + std::pow(h_y, 2))); // u_{i, j-1}
            }
            
            if (cur_y_id + 1 == y.size() - 1) {
                // u_{i, j+1} is already known from border conditions -> move the result to b
                b[cur_eq_id] += psi_1 * std::pow(h_x, 2) / (2 * (std::pow(h_x, 2) + std::pow(h_y, 2)));
            } else {
                A[cur_eq_id][mapping[cur_x_id][cur_y_id + 1]] = -std::pow(h_x, 2) / (2 * (std::pow(h_x, 2) + std::pow(h_y, 2))); // u_{i, j+1}
            }
            
            if (cur_x_id - 1 == 0) {
                // u_{i-1, j} is already known from border conditions -> move the result to b
                b[cur_eq_id] += phi_0 * std::pow(h_y, 2) / (2 * (std::pow(h_x, 2) + std::pow(h_y, 2)));
            } else {
                A[cur_eq_id][mapping[cur_x_id - 1][cur_y_id]] = -std::pow(h_y, 2) / (2 * (std::pow(h_x, 2) + std::pow(h_y, 2))); // u_{i-1, j}
            }
            
            if (cur_x_id + 1 == x.size() - 1) {
                // u_{i+1, j} is already known from border conditions -> move the result to b
                b[cur_eq_id] += phi_1 * std::pow(h_y, 2) / (2 * (std::pow(h_x, 2) + std::pow(h_y, 2)));
            } else {

                A[cur_eq_id][mapping[cur_x_id + 1][cur_y_id]] = -std::pow(h_y, 2) / (2 * (std::pow(h_x, 2) + std::pow(h_y, 2))); // u_{i+1, j}
            }
        }
    }
    
    // Step 3. Solve system of equations
    if (metod == 'iterative'){
	std::vector<std::vector<double>> ans = iterative(A, b, eps);
    
    	for (int cur_x_id = 1; cur_x_id < x.size() - 1; cur_x_id++) {
        	for (int cur_y_id = 1; cur_y_id < y.size() - 1; cur_y_id++) {
            	res[cur_x_id][cur_y_id] = ans[mapping[cur_x_id][cur_y_id]][0];
        }
     }
     if(metod == 'seidel'){
	std::vector<std::vector<double>> ans = seidel(A, b, eps);
    
    	for (int cur_x_id = 1; cur_x_id < x.size() - 1; cur_x_id++) {
        	for (int cur_y_id = 1; cur_y_id < y.size() - 1; cur_y_id++) {
            	res[cur_x_id][cur_y_id] = ans[mapping[cur_x_id][cur_y_id]][0];
        }
     }
     if (metod == 'relaxation'){
	std::vector<std::vector<double>> ans = relaxation(A, b, eps);
    
    	for (int cur_x_id = 1; cur_x_id < x.size() - 1; cur_x_id++) {
        	for (int cur_y_id = 1; cur_y_id < y.size() - 1; cur_y_id++) {
            	res[cur_x_id][cur_y_id] = ans[mapping[cur_x_id][cur_y_id]][0];
        }
     } 
    }
    
    return res;
}

std::vector<double> iterative(const std::vector<std::vector<double>>& A, const std::vector<double>& b, double eps) {
    int n = A.size();
    
    // Step 1. Ax = b -> x = alpha * x + beta
    std::vector<std::vector<double>> alpha(n, std::vector<double>(n, 0.0));
    std::vector<double> beta(n, 0.0);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                alpha[i][j] = 0.0;
            } else {
                alpha[i][j] = -A[i][j] / A[i][i];
            }
        }
        
        beta[i] = b[i] / A[i][i];
    }
    
    // Step 2. Iterating
    int iterations = 0;
    std::vector<double> cur_x(beta);
    bool converge = false;
    
    while (!converge) {
        std::vector<double> prev_x(cur_x);
        cur_x = std::vector<double>(n, 0.0);
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cur_x[i] += alpha[i][j] * prev_x[j];
            }
            
            cur_x[i] += beta[i];
        }
        
        double error = 0.0;
        
        for (int i = 0; i < n; i++) {
            error += std::abs(prev_x[i] - cur_x[i]);
        }
        
        iterations++;
        
        if (error <= eps) {
            converge = true;
        }
    }
    
    return cur_x;
}

std::vector<double> seidel_multiplication(const std::vector<std::vector<double>>& alpha, const std::vector<double>& x, const std::vector<double>& beta) {
    std::vector<double> res = x;
    
    for (int i = 0; i < alpha.size(); i++) {
        res[i] = beta[i];
        
        for (int j = 0; j < alpha[i].size(); j++) {
            res[i] += alpha[i][j] * res[j];
        }
    }
    
    return res;
}

std::vector<double> seidel(const std::vector<std::vector<double>>& A, const std::vector<double>& b, double eps) {
    int n = A.size();
    
    // Step 1. Ax = b -> x = alpha * x + beta
    std::vector<std::vector<double>> alpha(n, std::vector<double>(n, 0.0));
    std::vector<double> beta(n, 0.0);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                alpha[i][j] = 0.0;
            } else {
                alpha[i][j] = -A[i][j] / A[i][i];
            }
        }
        
        beta[i] = b[i] / A[i][i];
    }
    
    // Step 2. Iterating
    int iterations = 0;
    std::vector<double> cur_x = beta;
    bool converge = false;
    
    while (!converge) {
        std::vector<double> prev_x = cur_x;
        cur_x = seidel_multiplication(alpha, prev_x, beta);
        double error = 0.0;
        
        for (int i = 0; i < n; i++) {
            error += std::abs(prev_x[i] - cur_x[i]);
        }
        
        iterations++;
        
        if (error <= eps) {
            converge = true;
        }
    }
    
    return cur_x;
}


std::vector<double> relaxation(const std::vector<std::vector<double>>& A, const std::vector<double>& b, double eps, double w = 1.5) {
    int n = A.size();
    
    // Step 1. Ax = b -> x = alpha * x + beta
    std::vector<std::vector<double>> alpha(n, std::vector<double>(n, 0.0));
    std::vector<double> beta(n, 0.0);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                alpha[i][j] = 0.0;
            } else {
                alpha[i][j] = -A[i][j] / A[i][i];
            }
        }
        
        beta[i] = b[i] / A[i][i];
    }
    
    // Step 2. Iterating
    int iterations = 0;
    std::vector<double> cur_x = beta;
    bool converge = false;
    
    while (!converge) {
        std::vector<double> prev_x = cur_x;
        cur_x = seidel_multiplication(alpha, prev_x, beta);
        cur_x = w * cur_x + (1 - w) * prev_x;
        double error = 0.0;
        
        for (int i = 0; i < n; i++) {
            error += std::abs(prev_x[i] - cur_x[i]);
        }
        
        iterations++;
        
        if (error <= eps) {
            converge = true;
        }
    }
    
    return cur_x;
}

int main() {
    // Пример использования функции
    double x_begin = 0.0;
    double x_end = 1.05;
    double y_begin = 0.0;
    double y_end = 1.05;
    double h_x = 0.05;
    double h_y = 0.05;
    

    std::pair<double> x_range = std::make_pair(x_begin, x_end);
    std::pair<double, double> y_range = std::make_pair(y_begin, y_end);
    std::vector<std::vector<double>> analytical_solution = get_analytical_solution(x_range, y_range, h_x, h_y);
    
    // Вывод результатов
    for (int i = 0; i < analytical_solution.size(); i++) {
        for (int j = 0; j < analytical_solution[i].size(); j++) {
            std::cout << analytical_solution[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::vector<std::vector<double>>  iterative_solution;
    int relaxation_iters;
    std::tie(relaxation_solution,  iterative_iters) = finite_difference_schema(x_range, y_range, h_x, h_y,  iterative, 1e-7);
    solutions["relaxation solution"] =  iterative_solution;
    std::cout << "max abs error = " << max_abs_error( iterative_solution, analytical_solution) << std::endl;
    std::cout << "mean abs error = " << mean_abs_error( iterative_solution, analytical_solution) << std::endl;
    std::cout << "iterations = " <<  iterative_iters << std::endl;

    std::vector<std::vector<double>> seldel_solution;
    int relaxation_iters;
    std::tie(relaxation_solution, seldel_iters) = finite_difference_schema(x_range, y_range, h_x, h_y, seldel, 1e-7);
    solutions["relaxation solution"] = seldel_solution;
    std::cout << "max abs error = " << max_abs_error(seldel_solution, analytical_solution) << std::endl;
    std::cout << "mean abs error = " << mean_abs_error(seldel_solution, analytical_solution) << std::endl;
    std::cout << "iterations = " << seldel_iters << std::endl;


    std::vector<std::vector<double>> relaxation_solution;
    int relaxation_iters;
    std::tie(relaxation_solution, relaxation_iters) = finite_difference_schema(x_range, y_range, h_x, h_y, relaxation, 1e-7);
    solutions["relaxation solution"] = relaxation_solution;
    std::cout << "max abs error = " << max_abs_error(relaxation_solution, analytical_solution) << std::endl;
    std::cout << "mean abs error = " << mean_abs_error(relaxation_solution, analytical_solution) << std::endl;
    std::cout << "iterations = " << relaxation_iters << std::endl;

    return 0;
}
