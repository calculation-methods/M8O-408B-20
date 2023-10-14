#include <iostream>
#include <cmath>
#include <matplot/matplot.h>

using namespace std;

double l = M_PI / 2;

double psi(double x) {
    return 0;
}

double f(double x, double t) {
    return cos(x) * (cos(t) + sin(t));
}

double phi0(double t) {
    return sin(t);
}

double phi1(double t) {
    return -sin(t);
}

double solve(double x, double t) {
    return sin(t) * cos(x);
}

vector<double> solve_tri_diagonal_matrix(const vector<double> &lower_diagonal,
                                          const vector<double> &main_diagonal,
                                          const vector<double> &upper_diagonal,
                                          const vector<double> &free_coeffs) {
    int size = main_diagonal.size();
    vector<double> coefficientsP(size), coefficientsQ(size), solutionX(size);

    coefficientsP[0] = -upper_diagonal[0] / main_diagonal[0];
    coefficientsQ[0] = free_coeffs[0] / main_diagonal[0];

    for (int i = 1; i < size; i++) {
        double pTmp = -upper_diagonal[i] / (main_diagonal[i] + lower_diagonal[i] * coefficientsP[i - 1]);
        double qTmp = (free_coeffs[i] - lower_diagonal[i] * coefficientsQ[i - 1]) /
                      (main_diagonal[i] + lower_diagonal[i] * coefficientsP[i - 1]);
        coefficientsP[i] = pTmp;
        coefficientsQ[i] = qTmp;
    }

    solutionX[size - 1] = coefficientsQ[size - 1];
    for (int i = size - 2; i >= 0; i--) {
        solutionX[i] = coefficientsP[i] * solutionX[i + 1] + coefficientsQ[i];
    }

    return solutionX;
}

class task1 {
public:
    double alpha, beta, gamma, delta;
    string bound_type;
    double a, b, c, h, tau, sigma;

    task1(string bound_type, int N, int K, double T)
            : alpha(0), beta(1), gamma(1), delta(0), bound_type(bound_type), a(1), b(0), c(0) {
        h = l / N;
        tau = T / K;
        sigma = a * a * tau / (h * h);
    }

    vector<vector<double>> analytic(int N, int K, double T) {
        h = l / N;
        tau = T / K;
        vector<vector<double>> u(K, vector<double>(N, 0.0));

        for (int i = 0; i < K; i++) {
            for (int j = 0; j < N; j++) {
                u[i][j] = solve(j * h, i * tau);
            }
        }
        return u;
    }

    void set_coeffs_for_matrix_solution(vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d,
                   vector<vector<double>> &u, int k, int N, double T, int K) const {
        vector<double> t(K);
        for (int i = 0; i < K; i++) {
            t[i] = i * T / K;
        }

        for (int j = 1; j < N; j++) {
            a[j] = sigma;
            b[j] = -(1 + 2 * sigma);
            c[j] = sigma;
            d[j] = -u[k][j];
        }

        if (bound_type == "appr1") {
            a[0] = 0;
            b[0] = beta - (alpha / h);
            c[0] = alpha / h;
            d[0] = phi0(t[k]) / (beta - alpha / h);
            a[N - 1] = -gamma / h;
            b[N - 1] = gamma / h + delta;
            c[N - 1] = 0;
            d[N - 1] = phi1(t[k]) / (gamma / h + delta);
        } else if (bound_type == "appr2") {
            a[0] = 0;
            b[0] = -(1 + 2 * sigma);
            c[0] = sigma;
            d[0] = -(u[k - 1][0] + sigma * phi0(k * tau)) - tau * f(0, k * tau);
            a[N - 1] = sigma;
            b[N - 1] = -(1 + 2 * sigma);
            c[N - 1] = 0;
            d[N - 1] = -(u[k - 1][N - 1] + sigma * phi1(k * tau)) - tau * f((N - 1) * h, k * tau);
        } else if (bound_type == "appr3") {
            a[0] = 0;
            b[0] = -(1 + 2 * sigma);
            c[0] = sigma;
            d[0] = -((1 - sigma) * u[k - 1][1] + sigma / 2 * u[k - 1][0]) - tau * f(0, k * tau) -
                   sigma * phi0(k * tau);
            a[N - 1] = sigma;
            b[N - 1] = -(1 + 2 * sigma);
            c[N - 1] = 0;
            d[N - 1] = phi1(k * tau) + f((N - 1) * h, k * tau) * h / (2 * tau) * u[k - 1][N - 1];
        }
    }

    vector<vector<double>> implicit(int N, int K, double T) {
        vector<double> a(N, 0.0);
        vector<double> b(N, 0.0);
        vector<double> c(N, 0.0);
        vector<double> d(N, 0.0);
        vector<vector<double>> u(K, vector<double>(N, 0.0));

        for (int i = 1; i < N - 1; i++) {
            u[0][i] = psi(i * h);
        }
        u[0][N - 1] = 0;

        for (int k = 1; k < K; k++) {
            set_coeffs_for_matrix_solution(a, b, c, d, u, k, N, T, K);
            u[k] = solve_tri_diagonal_matrix(a, b, c, d);
        }

        return u;
    }

    vector<vector<double>> eexplicit(int N, int K, double T) const {
        vector<vector<double>> u(K, vector<double>(N, 0.0));
        vector<double> t(K);
        vector<double> x(N);

        for (int i = 0; i < K; i++) {
            t[i] = i * T / K;
        }

        for (int i = 0; i < N; i++) {
            x[i] = i * M_PI / 2 / N;
        }

        for (int j = 1; j < N - 1; j++) {
            u[0][j] = psi(j * h);
        }

        for (int k = 1; k < K; k++) {
            for (int j = 1; j < N - 1; j++) {
                u[k][j] = (u[k - 1][j + 1] * (a * a * tau / (h * h))
                           - 2 * u[k - 1][j] * (a * a * tau / (h * h))
                           + u[k - 1][j - 1] * (a * a * tau / (h * h))
                           + u[k - 1][j]
                           + tau * f(x[j], t[k]));
            }

            if (bound_type == "appr1") {
                u[k][0] = phi0(t[k]);
                u[k][N - 1] = (phi1(t[k]) + gamma / h * u[k][N - 2]) / (delta + gamma / h);
            } else if (bound_type == "appr2") {
                u[k][0] = phi0(t[k]);
                u[k][N - 1] = (((2.0 * gamma * a / h / (2.0 * a + h * b) * u[k][N - 2]) +
                                (gamma * h / tau / (2.0 * a + h * b) * u[k - 1][N - 1]) +
                                (gamma * h * c / (2.0 * a + h * b) * f(l, t[k]) + phi1(t[k])) /
                                (2.0 * gamma * a / h / (2.0 * a + h * b) + gamma * h / tau / (2.0 * a + h * b) -
                                 gamma * h * c / (2.0 * a + h * b) * c + delta)));
            } else if (bound_type == "appr3") {
                u[k][0] = phi0(t[k]);
                u[k][N - 1] = (phi1(t[k]) + u[k][N - 2] / h + 2 * tau * u[k - 1][N - 1] / h) /
                              (1 / h + 2 * tau / h);
            }
        }

        return u;
    }

    vector<vector<double>> niklson(int N, int K, double T) {
        double theta = 0.5;
        vector<vector<double>> u(K, vector<double>(N, 0.0));
        vector<double> a(N, 0.0);
        vector<double> b(N, 0.0);
        vector<double> c(N, 0.0);
        vector<double> d(N, 0.0);

        for (int i = 1; i < N - 1; i++) {
            u[0][i] = psi(i * h);
        }

        for (int k = 1; k < K; k++) {
            set_coeffs_for_matrix_solution(a, b, c, d, u, k, N, T, K);

            vector<double> tmp_imp = solve_tri_diagonal_matrix(a, b, c, d);
            vector<double> tmp_exp(N, 0.0);

            tmp_exp[0] = phi0(k * tau);
            for (int j = 1; j < N - 1; j++) {
                tmp_exp[j] = sigma * u[k - 1][j + 1] + (1 - 2 * sigma) * u[k - 1][j] +
                             sigma * u[k - 1][j - 1] + tau * f(j * h, k * tau);
            }
            tmp_exp[N - 1] = phi1(k * tau);

            for (int j = 0; j < N; j++) {
                u[k][j] = theta * tmp_imp[j] + (1 - theta) * tmp_exp[j];
            }
        }

        return u;
    }


};

vector<vector<double>> compare_error(const vector<vector<double>> &numerical, const vector<vector<double>> &analytic) {
    vector<vector<double>> error;

    for (int i = 0; i < numerical.size(); i++) {
        vector<double> row;
        for (int j = 0; j < numerical[i].size(); j++) {
            row.push_back(abs(numerical[i][j] - analytic[i][j]));
        }
        error.push_back(row);
    }

    return error;
}

void plot_results(map<string, vector<vector<double>>> &dict, int time) {
    using namespace matplot;

    tiledlayout(2, 2);
    nexttile();
//    ylim({-0.3, 0.3});
    title("Линии уровня");
    auto p = plot(dict["implicit"][time]);
    p->display_name("implicit");
    hold(on);
    p = plot(dict["explicit"][time]);
    p->display_name("explicit");
    p = plot(dict["niklson"][time]);
    p->display_name("niklson");
    p = plot(dict["analytic"][time]);
    p->display_name("analytic");
    matplot::legend();
    hold(off);


    nexttile();
    title("Погрешность explicit");
    plot(compare_error(dict["explicit"], dict["analytic"])[time]);
    ylabel("Err");
    xlabel("t");

    nexttile();
    title("Погрешность implicit");
    plot(compare_error(dict["implicit"], dict["analytic"])[time]);
    ylabel("Err");
    xlabel("t");

    nexttile();
    title("Погрешность niklson");
    plot(compare_error(dict["niklson"], dict["analytic"])[time]);
    ylabel("Err");
    xlabel("t");
}

int main() {
    int N = 12;
    int K = 10000;
    int T = 1;

    for (string bound_type: {"appr1", "appr2", "appr3"}) {
        auto task = task1(bound_type, N, K, T);

        map<string, vector<vector<double>>> results;

        results["implicit"] = task.implicit(N, K, T);
        results["explicit"] = task.eexplicit(N, K, T);
        results["niklson"] = task.niklson(N, K, T);
        results["analytic"] = task.analytic(N, K, T);
        matplot::figure();
        plot_results(results, 2);
        matplot::show();
    }

    return 0;
}
