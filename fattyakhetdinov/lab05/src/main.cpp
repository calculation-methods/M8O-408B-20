#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using std::ofstream;

const int N = 10;
const int K = 100;
const int T = 1;
const double L = M_PI / 2;
const double H = L / N;
const double TAU = static_cast<double>(T) / static_cast<double>(K);
const double SIGMA = TAU / pow(H, 2);

template <typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

    std::vector<double> linspace_vector;
    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double amount = static_cast<double>(num_in);
    if (amount == 0)
    {
        return linspace_vector;
    }
    else if (amount == 1)
    {
        linspace_vector.push_back(start);
        return linspace_vector;
    }
    double delta = (end - start) / (amount - 1);
    for (int i = 0; i < amount - 1; ++i)
    {
        linspace_vector.push_back(start + delta * i);
    }
    linspace_vector.push_back(end);
    return linspace_vector;
}

void debug_print_vector(std::vector<double> vec, std::string name)
{
    std::cout << "------------------" << std::endl
              << "[DEBUG] " << name << " size: " << vec.size() << std::endl;
    for (double d : vec)
        std::cout << d << " ";
    std::cout << std::endl
              << "------------------" << std::endl;
}
void debug_print_vector2(std::vector<std::vector<double>> vec, std::string name)
{
    std::cout << "------------------" << std::endl
              << "[DEBUG] " << name << " size: " << vec.size() << "x" << vec[0].size() << std::endl;
    for (auto v : vec)
    {
        for (auto elem : v)
        {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    std::cout
        << "------------------" << std::endl;
}
void save_vector2_to_file(std::vector<std::vector<double>> vec, std::string filename)
{
    ofstream file(filename);
    for (auto v : vec)
    {
        for (auto elem : v)
        {
            file << elem << " ";
        }
        file << "\n";
    }
    std::cout << "Saved to: " << filename << std::endl;
}
template <typename T>
void debug_print(T value, std::string name)
{
    std::cout << "[DEBUG] " << name << ": " << value << std::endl;
}

double PSI(double x)
{
    return 0;
}

double F(double x, double t)
{
    return cos(x) * (cos(t) + sin(t));
}

double PHI0(double t)
{
    return sin(t);
}
double PHI1(double t)
{
    return -sin(t);
}
double SOLUTION(double x, double t)
{
    return sin(t) * cos(x);
}

std::vector<std::vector<double>> analitical_solution(std::vector<double> x, std::vector<double> t)
{
    std::vector<std::vector<double>> res;
    for (int i = 0; i < K; ++i)
    {
        std::vector<double> subres(N);
        res.push_back(subres);
    }
    for (int i = 0; i < K; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            res[i][j] = SOLUTION(x[j], t[i]);
        }
    }
    return res;
}

std::vector<std::vector<double>> explicit_solution(int app)
{
    std::vector<std::vector<double>> res;
    for (int i = 0; i < K; ++i)
    {
        std::vector<double> subres(N);
        res.push_back(subres);
    }
    for (int i = 0; i < N; ++i)
    {
        res[0][i] = PSI(i * H);
    }
    for (int i = 0; i < K; ++i)
    {
        res[i][0] = PHI0(i * TAU);
    }
    for (int i = 1; i < K; ++i)
    {
        for (int j = 1; j < N - 1; ++j)
        {
            res[i][j] = res[i - 1][j] + TAU * (res[i - 1][j - 1] - 2 * res[i - 1][j] + res[i - 1][j + 1]) / pow(H, 2) + TAU * F(j * H, TAU * i);
        }
    }
    if (app == 1)
    {
        for (int i = 0; i < K; ++i)
        {
            res[i][res[i].size() - 1] = PHI1(TAU * i) * H + res[i][res[i].size() - 2];
        }
    }
    else if (app == 2)
    {
        for (int i = 0; i < K; ++i)
        {
            res[i][res[i].size() - 1] = (PHI1(i * TAU) * 2 * H - res[i][res[i].size() - 3] + 4 * res[i][res[i].size() - 2]) / 3;
        }
    }
    else if (app == 3)
    {
        for (int i = 1; i < K; ++i)
        {
            res[i][res[i].size() - 1] = (PHI1(i * TAU) + res[i][res[i].size() - 2] / H + 2 * TAU * res[i - 1][res[i].size() - 1] / H) / (1 / H + 2 * TAU / H);
        }
    }
    return res;
}

std::vector<double> solve_tridiagonal(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d)
{
    std::vector<double> p, q;
    std::vector<double> res(a.size());
    p.push_back(-c[0] / b[0]);
    q.push_back(d[0] / b[0]);
    for (int i = 1; i < a.size(); ++i)
    {
        p.push_back(-c[i] / (b[i] + a[i] * p[i - 1]));
        q.push_back((d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]));
    }

    res[res.size() - 1] = q[q.size() - 1];
    for (int i = res.size() - 2; i > -1; --i)
    {
        res[i] = p[i] * res[i + 1] + q[i];
    }

    return res;
}

std::vector<std::vector<double>> implicit_solution()
{
    std::vector<double> a(N);
    std::vector<double> b(N);
    std::vector<double> c(N);
    std::vector<double> d(N);
    std::vector<std::vector<double>> res;
    for (int i = 0; i < K; ++i)
    {
        std::vector<double> subres(N);
        res.push_back(subres);
    }
    for (int i = 0; i < N; ++i)
    {
        res[0][i] = PSI(i * H);
    }
    for (int i = 0; i < K; ++i)
    {
        res[i][0] = PHI0(i * TAU);
    }
    for (int i = 1; i < K; ++i)
    {
        a[0] = 0;
        b[0] = -(1 + 2 * SIGMA);
        c[0] = SIGMA;
        d[0] = -res[i - 1][0] - SIGMA * PHI0(i * TAU);
        for (int j = 1; j < N - 1; ++j)
        {
            a[j] = SIGMA;
            b[j] = -(1 + 2 * SIGMA);
            c[j] = SIGMA;
            d[j] = -res[i - 1][j] - TAU * F(j * H, (i - 1) * TAU);
        }
        a[a.size() - 1] = SIGMA;
        b[b.size() - 1] = -(1 + 2 * SIGMA);
        c[c.size() - 1] = 0;
        d[d.size() - 1] = -H * PHI1(TAU * i) * H - res[i][res[i].size() - 1] - TAU * F(N * H, TAU * (i + 1));
        // debug_print_vector(a, "a");
        // debug_print_vector(b, "b");
        // debug_print_vector(c, "c");
        // debug_print_vector(d, "d");
        // exit(1);
        res[i] = solve_tridiagonal(a, b, c, d);
    }
    return res;
}

std::vector<std::vector<double>> crank_nicholson_solution()
{
    double tetta = 0.5;
    std::vector<std::vector<double>> res;
    for (int i = 0; i < K; ++i)
    {
        std::vector<double> subres(N);
        res.push_back(subres);
    }
    std::vector<std::vector<double>> implicit_sol = implicit_solution();
    std::vector<std::vector<double>> explicit_sol = explicit_solution(1);
    // debug_print_vector2(implicit_sol, "implicit");
    // debug_print_vector2(explicit_sol, "explicit");
    // debug_print_vector2(res, "res");
    for (int i = 0; i < K; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            res[i][j] = implicit_sol[i][j] * tetta + explicit_sol[i][j] * (1 - tetta);
        }
    }
    return res;
}

int main()
{
    char debug;
    std::cout << "Debug mode (all vectors will be printed to stdout)? y/n" << std::endl;
    std::cin >> debug;
    std::vector<double> x = linspace(0., L, N);
    std::vector<double> t = linspace(0, T, K);

    auto analytical_sol = analitical_solution(x, t);
    auto explicit_sol_1 = explicit_solution(1);
    auto explicit_sol_2 = explicit_solution(2);
    auto explicit_sol_3 = explicit_solution(3);
    auto implicit_sol = implicit_solution();
    auto crank_nicholson_sol = crank_nicholson_solution();
    save_vector2_to_file(analytical_sol, "analytical_sol.txt");
    save_vector2_to_file(explicit_sol_1, "explicit_sol_1.txt");
    save_vector2_to_file(explicit_sol_2, "explicit_sol_2.txt");
    save_vector2_to_file(explicit_sol_3, "explicit_sol_3.txt");
    save_vector2_to_file(implicit_sol, "implicit_sol.txt");
    save_vector2_to_file(crank_nicholson_sol, "crank_nicholson_sol.txt");

    if (debug == 'y')
    {
        debug_print_vector(x, "x linspace");
        debug_print_vector(t, "t linspace");
        debug_print(L, "l");
        debug_print(SIGMA, "sigma");
        debug_print(TAU, "tau");
        debug_print(H, "h");
        debug_print_vector2(analytical_sol, "Analytical solution");
        debug_print_vector2(explicit_sol_1, "Explicit one point");
        debug_print_vector2(explicit_sol_2, "Explicit two points");
        debug_print_vector2(explicit_sol_3, "Explicit three points");
        debug_print_vector2(implicit_sol, "Implicit solution");
        debug_print_vector2(crank_nicholson_sol, "Crank-Nicholson solution");
    }
}
