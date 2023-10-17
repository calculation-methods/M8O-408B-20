#include <iostream>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <stdint.h>

using namespace std;


const double PI = acos(0.) * 2;

double f(double x, double t) {
    return 0;
}

double alpha = 1;
double beta = 1;

double g0(double a, double b, double c, double t) {
    return exp((c - a) * t) * (cos(b * t) + sin(b * t));
}

double g1(double a, double b, double c,  double t) {
    return -exp((c - a) * t)* (cos(b * t) + sin(b * t));
}

vector<double> ThreeDiagonal(int n, vector<double> m, vector<double> d) {
    vector<double> x(n);
    vector<double> q(n);
    vector<double> p(n);

    p[0] = -m[1] / m[0];
    q[0] = d[0] / m[0];
    
    int k = 1;
    for (uint i = 2; i < m.size() - 2; i += 3) {
        p[k] = - m[i + 2] / (m[i + 1] + m[i] * p[k - 1]);
        q[k] = (d[k] - m[i] * q[k - 1])/(m[i + 1] + m[i] * p[k - 1]);
        k += 1;
    }

    p[k] = 0;
    q[k] = (d[k] - m[m.size() - 2] * q[k - 1]) / (m[m.size() - 1] + m[m.size() - 2] * p[k - 1]);

    x[n - 1] = q[n - 1];

    for (int i = n - 2; i >= 0; --i) {
        x[i] = p[i] * x[i + 1] + q[i];
    }

    return x;
}

vector<vector<double>> ThreeDiagonalsolutuionsving(vector<vector<double>> u, int n, int k, double a, double c, double b, double tau, double h, int t = 0) {
    

    for (int j = 1; j < k; ++j) {

        vector<double> m((n - 2) * 3 + 4);
        vector<double> d(n);

        if (t == 0) {
            m[0] = -alpha / h + beta;
            m[1] = alpha / h;

            d[0] = g0(a, b, c, tau * (j + 1));
            d[n - 1] = g1(a, b, c, tau * (j + 1));

            m[m.size() - 2] = -alpha / h;
            m[m.size() - 1] = alpha / h + beta;
        }

        int mi = 2;


        for (int i = 1; i < n - 1; ++i) {
            m[mi] = (a * tau / (h * h)  - b * tau / (2 * h));
            m[mi + 1] = -1 - 2 * a * tau / (h * h) + c * tau;
            m[mi + 2] = (a * tau / (h * h)  + b * tau / (2 * h));
            mi += 3;
            d[i] = -u[i][j - 1] - tau * f(h * i, j * tau);
            
        }

        if (t == 1) {
            double c = (- alpha / (2 * h)) / m[4];
            m[0] = (-3 * alpha / (2 * h)) + beta - c * m[2];
            m[1] = alpha * 2 / h - c * m[3];
            d[0] = g0(a, b, c, (j + 1) * tau) - c * d[1];

            c = (alpha / (2 * h))  / m[m.size() - 5];
            m[m.size() - 2] = -alpha * 2 / h - c * m[m.size() - 4];
            m[m.size() - 1] = (3 * alpha / (2 * h)) + beta - m[m.size() - 3] * c;
            d[n - 1] = g1(a, b, c, (j + 1) * tau) -  d[n - 2] * c;

        }

        if (t == 2) {

            double c = h - b * h * h / (2 * a);
            m[0] = -(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / c * alpha + beta;
            m[1] = alpha / c;

            d[0] = g0(a, b, c, (j + 1) * tau) - (u[0][j - 1] * alpha * h * h / (2 * a * tau) + alpha * f(0, j * tau) * h * h / (2 * a)) / c;


            c = -h - b * h * h / (2 * a);
            m[m.size() - 2] = alpha / c;
            m[m.size() - 1] = -(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / c * alpha + beta;

            d[n - 1] = g1(a, b, c, (j + 1) * tau) - (u[n - 1][j - 1] * alpha * h * h / (2 * a * tau) + alpha * f(0, j * tau) * h * h / (2 * a)) / c;
        }

        d = ThreeDiagonal(n, m, d);

        for (int i = 0; i < n; ++i) {
            u[i][j] = d[i];
        }
    }

    return u;

}

vector<vector<double>> MethodOfFiniteDifference(vector<vector<double>> u, int n, int k, double a, double c, double b, double tau, double h, int t) {

    for (int j = 0; j < k - 1; ++j) {
        for (int i = 1; i < n - 1; ++i) {
            u[i][j + 1] = tau * (a * (u[i - 1][j] - 2 * u[i][j] + u[i + 1][j]) / (h * h) + 
                          b * (u[i + 1][j] - u[i - 1][j]) / (2 * h) + 
                          c * u[i][j] + f(h * i, j * tau)) + u[i][j];

            if (t == 0) {
                u[0][j + 1] = (g0(a, b, c, (j + 1) * tau) - alpha / h * u[1][j + 1]) / (-alpha / h + beta);
                u[n - 1][j + 1] = (g1(a, b, c, (j + 1) * tau) + alpha / h * u[n - 2][j + 1]) / (alpha / h + beta);
            }

            if (t == 1) {
                u[0][j + 1] = (g0(a, b, c, (j + 1) * tau) - alpha * 2 / h * u[1][j + 1] + alpha / (2 * h) * u[2][j + 1]) / (-3 * alpha / ( 2 * h) + beta);
                u[n - 1][j + 1] = (g1(a, b, c, (j + 1) * tau) - alpha / (h * 2) * u[n - 3][j + 1] + 2 * alpha / h * u[n - 2][j + 1]) / (3 * alpha / (h * 2) + beta);
            }

            if (t == 2) {
                double c = h - b * h * h / (2 * a);
                u[0][j + 1] = (g0(a, b, c, (j + 1) * tau) - alpha / c * u[1][j + 1] - (alpha * h * h / (2 * tau * a) * u[0][j] + alpha *  f(0, j * tau) * h * h / (2 * a)) / c) / 
                              (-(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / c * alpha + beta);

                c = -h - b * h * h / (2 * a);
                u[n - 1][j + 1] = (g1(a, b, c, (j + 1) * tau) - alpha / c * u[n - 2][j + 1] - (alpha * h * h / (2 * tau * a) * u[n - 1][j] + alpha *  f(0, j * tau) * h * h / (2 * a)) / c) / 
                              (-(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / c * alpha + beta);
            }
        }
    }

    return u;

}

vector<vector<double>> MainsolutuionsvingMethod(vector<vector<double>> u, int n, int k, double a, double c, double b, double tau, double h, double theta, int t) {
    for (int j = 1; j < k; ++j) {

        vector<double> m((n - 2) * 3 + 4);
        vector<double> d(n);

        if (t == 0) {
            m[0] = -alpha / h + beta;
            m[1] = alpha / h;

            d[0] = g0(a, b, c, tau * (j + 1));
            d[n - 1] = g1(a, b, c, tau * (j + 1));

            m[m.size() - 2] = -alpha / h;
            m[m.size() - 1] = alpha / h + beta;
        }

        int mi = 2;


        for (int i = 1; i < n - 1; ++i) {
            m[mi] = theta * (a * tau / (h * h)  - b * tau / (2 * h));
            m[mi + 1] = -1 + theta * (- 2 * a * tau / (h * h) + c * tau);
            m[mi + 2] = theta * (a * tau / (h * h)  + b * tau / (2 * h));
            mi += 3;


            d[i] = -u[i][j - 1]  - tau * (theta * f(h * i, j * tau) + (1 - theta) * ((u[i - 1][j - 1] - 2 * u[i][j - 1] + u[i + 1][j - 1]) / (h * h) * a  + 
                (u[i + 1][j - 1] - u[i - 1][j - 1]) / (2 * h) * b + 
                c * u[i][j - 1]));
        }

          if (t == 1) {
            double c = (- alpha / (2 * h)) / m[4];
            m[0] = (-3 * alpha / (2 * h)) + beta - c * m[2];
            m[1] = alpha * 2 / h - c * m[3];
            d[0] = g0(a, b, c, (j + 1) * tau) - c * d[1];

            c = (alpha / (2 * h))  / m[m.size() - 5];
            m[m.size() - 2] = -alpha * 2 / h - c * m[m.size() - 4];
            m[m.size() - 1] = (3 * alpha / (2 * h)) + beta - m[m.size() - 3] * c;
            d[n - 1] = g1(a, b, c, (j + 1) * tau) -  d[n - 2] * c;
        }


        if (t == 2) {

            double c = h - b * h * h / (2 * a);
            m[0] = -(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / c * alpha + beta;
            m[1] = alpha / c;

            d[0] = g0(a, b, c, (j + 1) * tau) - (u[0][j - 1] * alpha * h * h / (2 * a * tau) + alpha * f(0, j * tau) * h * h / (2 * a)) / c;


            c = -h - b * h * h / (2 * a);
            m[m.size() - 2] = alpha / c;
            m[m.size() - 1] = -(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / c * alpha + beta;

            d[n - 1] = g1(a, b, c, (j + 1) * tau) - (u[n - 1][j - 1] * alpha * h * h / (2 * a * tau) + alpha * f(0, j * tau) * h * h / (2 * a)) / c;
        }

        d = ThreeDiagonal(n, m, d);

        for (int i = 0; i < n; ++i) {
            u[i][j] = d[i];
        }

    }

    return u;

}


int main() {

    double a, b, c, sigma, l, T, tau, h, theta;

    int n, k;
    cin >> l >> n >> T >> sigma >> a >> b >> c >>theta;

    h = l / (n - 1);
    tau = sigma * h * h / a;
    k = T / tau + 1;


    vector<vector<double>> u(n, vector<double>(k));
    for (int i = 0; i < n; ++i) {
        u[i][0] = sin(i * h);
    }
    vector<vector<double>> solutuions(n, vector<double>(k));


    FILE * ResultFile;

    if ((ResultFile = fopen("data.txt", "w")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }

    fprintf(ResultFile,"%s", (to_string(a) + "\n").c_str());
    fprintf(ResultFile, "%s", (to_string(b) + "\n").c_str());
    fprintf(ResultFile, "%s", (to_string(c) + "\n").c_str());
    fprintf(ResultFile, "%s", (to_string(tau) + "\n").c_str());
    fprintf(ResultFile, "%s", (to_string(h) + "\n").c_str());
    fprintf(ResultFile, "%s", (to_string(k) + "\n").c_str());
    fprintf(ResultFile, "%s", (to_string(n) + "\n").c_str());


    solutuions = MethodOfFiniteDifference(u, n, k, a, c, b, tau, h, 0);


    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(ResultFile, "%s", (to_string(solutuions[i][j]) + "\n").c_str());
        }
    }


    solutuions = ThreeDiagonalsolutuionsving(u, n, k, a, c, b, tau, h, 0);



    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(ResultFile, "%s", (to_string(solutuions[i][j]) + "\n").c_str());
        }
    }

    solutuions = MainsolutuionsvingMethod(u, n, k, a, c, b, tau, h, theta, 0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(ResultFile, "%s", (to_string(solutuions[i][j]) + "\n").c_str());
        }
    }


    solutuions = MethodOfFiniteDifference(u, n, k, a, c, b, tau, h, 1);


    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(ResultFile, "%s", (to_string(solutuions[i][j]) + "\n").c_str());
        }
    }


    solutuions = ThreeDiagonalsolutuionsving(u, n, k, a, c, b, tau, h, 1);



    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(ResultFile, "%s", (to_string(solutuions[i][j]) + "\n").c_str());
        }
    }

    solutuions = MainsolutuionsvingMethod(u, n, k, a, c, b, tau, h, theta, 1);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(ResultFile, "%s", (to_string(solutuions[i][j]) + "\n").c_str());
        }
    }


    solutuions = MethodOfFiniteDifference(u, n, k, a, c, b, tau, h, 2);


    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(ResultFile, "%s", (to_string(solutuions[i][j]) + "\n").c_str());
        }
    }

    solutuions = ThreeDiagonalsolutuionsving(u, n, k, a, c, b, tau, h, 2);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(ResultFile, "%s", (to_string(solutuions[i][j]) + "\n").c_str());
        }
    }

    solutuions = MainsolutuionsvingMethod(u, n, k, a, c, b, tau, h, theta, 2);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(ResultFile, "%s", (to_string(solutuions[i][j]) + "\n").c_str());
        }
    }

    fclose(ResultFile);
}