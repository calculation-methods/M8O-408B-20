#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>

const double MAX_X = M_PI / 2;
const double MAX_Y = M_PI;
const double MAX_T = 1;

const double A = 1;
const double B = 1;
const double MU = 1;

// y = 0
double phi_1(double x, double t)
{
    return 0;
}

// y = max
double phi_d2(double x, double t)
{
    return -sin(x) * sin(MU * t);
}

// x = 0
double phi_3(double y, double t)
{
    return 0;
}

// x = max
double phi_4(double y, double t)
{
    return sin(y) * sin(MU * t);
}

// t = 0
double psi(double x, double y)
{
    return 0;
}

double f(double x, double y, double t)
{
    return sin(x) * sin(y) * (MU * cos(MU * t) + (A + B) * sin(MU * t));
}

double correct(double x, double y, double t)
{
    return sin(x) * sin(y) * sin(MU * t);
}

std::vector<double> progon(std::vector<double> ar, std::vector<double> br, std::vector<double> cr, std::vector<double> dr)
{
    int len = br.size();
    if (dr.size() != len)
    {
        throw std::invalid_argument("bad array dimensions");
    }
    if (ar.size() + 1 == len)
    {
        ar.insert(ar.begin(), 0);
    }
    else if (ar.size() != len)
    {
        throw std::invalid_argument("bad array dimensions");
    }
    if (cr.size() + 1 == len)
    {
        cr.push_back(0); // Corrected
    }
    else if (cr.size() != len)
    {
        throw std::invalid_argument("bad array dimensions");
    }
    std::vector<double> P(len, 0);
    std::vector<double> Q(len, 0);

    // Прямой ход
    P[0] = -cr[0] / br[0];
    Q[0] = dr[0] / br[0];
    for (int i = 1; i < len; ++i)
    {
        double denom = br[i] + ar[i] * P[i - 1];
        P[i] = -cr[i] / denom;
        Q[i] = (dr[i] - ar[i] * Q[i - 1]) / denom;
    }

    // Обратный ход
    std::vector<double> result(len);
    result[len - 1] = Q[len - 1];
    for (int i = len - 2; i >= 0; --i)
    {
        result[i] = P[i] * result[i + 1] + Q[i];
    }
    return result;
}

class griddy
{
    int max_x;
    int max_y;
    int max_ti;
    double step_x;
    double step_y;
    double step_ti;
    std::vector<double> vals;

    void process_arguments(int x, int y, int t) const
    {
        if ((x < 0) || (x >= max_x))
        {
            throw std::invalid_argument("индекс x вне диапазона");
        }
        if ((y < 0) || (y >= max_y))
        {
            throw std::invalid_argument("индекс y вне диапазона");
        }
        if ((t < 0) || (t > max_ti))
        {
            throw std::invalid_argument("индекс t вне диапазона");
        }
    }
    double get_x(int x) const
    {
        process_arguments(x, 0, 0);
        return step_x * x;
    }
    double get_y(int y) const
    {
        process_arguments(0, y, 0);
        return step_y * y;
    }
    double get_t(int ti) const
    {
        process_arguments(0, 0, ti);
        return step_ti * ti;
    }

    void halfstep_unclear_x(int j, int k)
    {
        std::vector<double> as(max_x - 2, 0);
        std::vector<double> bs(max_x - 2, 0);
        std::vector<double> cs(max_x - 2, 0);
        std::vector<double> ds(max_x - 2, 0);

        for (int i = 1; i < max_x - 1; i++)
        {
            double sa = -A / pow(step_x, 2);
            double sb = (1.0 / step_ti) + (2.0 * A / pow(step_x, 2));
            double sc = -A / pow(step_x, 2);
            double ff = f(get_x(i), get_y(j), get_t(k - 1));
            double sd = (ff / 2) + (U(i, j, k - 1) / step_ti);
            if (i == 1)
            {
                sd -= U(0, j, k) * sa;
                sa = 0;
            }
            if (i == (max_x - 2))
            {
                sd -= U(max_x - 1, j, k) * sc;
                sc = 0;
            }
            as[i - 1] = sa;
            bs[i - 1] = sb;
            cs[i - 1] = sc;
            ds[i - 1] = sd;
        }
        std::vector<double> rez = progon(as, bs, cs, ds);
        for (int i = 1; i < max_x - 1; i++)
        {
            U_mut(i, j, k) = rez[i - 1];
        }
    }

    void halfstep_unclear_y(int i, int k)
    {
        std::vector<double> as(max_y - 2, 0);
        std::vector<double> bs(max_y - 2, 0);
        std::vector<double> cs(max_y - 2, 0);
        std::vector<double> ds(max_y - 2, 0);

        for (int j = 1; j < max_y - 1; j++)
        {
            double sa = -B / pow(step_y, 2);
            double sb = (1.0 / step_ti) + (2.0 * B / pow(step_y, 2));
            double sc = -B / pow(step_y, 2);
            double ff = f(get_x(i), get_y(j), get_t(k - 1));
            double sd = (ff / 2) + (U(i, j, k - 1) / step_ti);
            if (j == 1)
            {
                sd -= U(i, 0, k) * sa;
                sa = 0;
            }
            if (j == (max_y - 2))
            {
                sd -= U(i, max_y - 1, k) * sc;
                sc = 0;
            }
            as[j - 1] = sa;
            bs[j - 1] = sb;
            cs[j - 1] = sc;
            ds[j - 1] = sd;
        }
        std::vector<double> rez = progon(as, bs, cs, ds);
        for (int j = 1; j < max_y - 1; j++)
        {
            U_mut(i, j, k) = rez[j - 1];
        }
    }

public:
    griddy(int x_steps, int y_steps, int t_steps)
    {
        max_x = x_steps;
        max_y = y_steps;
        max_ti = t_steps + 1;
        step_x = MAX_X / (max_x - 1);
        step_y = MAX_Y / (max_y - 1);
        step_ti = MAX_T / (max_ti - 1);
        vals = std::vector<double>(max_x * max_y * max_ti);
    }
    griddy(const griddy &other)
    {
        max_x = other.max_x;
        max_y = other.max_y;
        max_ti = other.max_ti;
        step_x = other.step_x;
        step_y = other.step_y;
        step_ti = other.step_ti;
        vals = std::vector<double>(other.vals.size());
        for (int i = 0; i < max_x * max_y * max_ti; i++)
        {
            vals[i] = other.vals[i];
        }
    }
    void annul()
    {
        for (int i = 0; i < vals.size(); i++)
        {
            vals[i] = 0;
        }
    }
    int idx(int i, int j, int k) const
    {
        return (k * max_ti + j) * max_y + i;
        // return (i * max_x + j) * max_y + k;
    }
    double &U_mut(int i, int j, int k)
    {
        process_arguments(i, j, k);
        return vals[idx(i, j, k)];
    }
    double U(int i, int j, int k) const
    {
        process_arguments(i, j, k);
        return vals[idx(i, j, k)];
    }
    void cheat_set_correct()
    {
        for (int i = 0; i < max_x; i++)
        {
            for (int j = 0; j < max_y; j++)
            {
                for (int k = 0; k < max_ti; k++)
                {
                    U_mut(i, j, k) = correct(get_x(i), get_y(j), get_t(k));
                }
            }
        }
    }
    void set_start()
    {
        for (int i = 0; i < max_x; i++)
        {
            for (int j = 0; j < max_y; j++)
            {
                U_mut(i, j, 0) = 0; // psi(get_x(i), get_y(j));
            }
        }
        for (int k = 0; k < max_ti; k++)
        {
            for (int i = 0; i < max_x; i++)
            {
                U_mut(i, 0, k) = phi_1(get_x(i), get_t(k));
            }
            for (int j = 0; j < max_y; j++)
            {
                U_mut(0, j, k) = phi_3(get_y(j), get_t(k));
                U_mut(max_x - 1, j, k) = phi_4(get_y(j), get_t(k));
            }
        }
    }
    double square_error()
    {
        double result = 0;
        for (int i = 0; i < vals.size(); i++)
        {
            result += pow(vals[i], 2);
        }
        return sqrt(result);
    }
    double n_row_error(int k)
    {
        double res = 0;
        for (int i = 0; i < max_x; i++)
        {
            for (int j = 0; j < max_y; j++)
            {
                res += pow(U(i, j, k), 2);
            }
        }
        return sqrt(res);
    }
    griddy diff(const griddy &other) const
    {
        if ((other.max_x != max_x) || (other.max_y != max_y) || (other.max_ti != max_ti))
        {
            throw std::invalid_argument("размеры не совпадают");
        }
        griddy rez(max_x, max_y, max_ti);
        for (int i = 0; i < max_x; i++)
        {
            for (int j = 0; j < max_y; j++)
            {
                for (int k = 0; k < max_ti; k++)
                {
                    double diff_value = this->U(i, j, k) - other.U(i, j, k);
                    *(rez.vals.begin() + (k * max_x * max_y + j * max_x + i)) = diff_value;
                }
            }
        }
        return rez;
    }

    void print_layer(int k)
    {
        for (int j = 0; j < max_y; j++)
        {
            for (int i = 0; i < max_x; i++)
            {
                printf("%5.2lf ", U(i, j, k));
            }
            printf("\n");
        }
    }
    void print_all()
    {
        for (int k = 0; k < max_ti; k++)
        {
            print_layer(k);
            printf("\n");
        }
        printf("-------------------------\n");
    }

    void full_halfstep_unclear_x(int k)
    {
        for (int j = 0; j < max_y; j++)
        {
            halfstep_unclear_x(j, k);
        }
    }

    void full_halfstep_unclear_y(int k)
    {
        for (int i = 0; i < max_x; i++)
        {
            halfstep_unclear_y(i, k);
        }
    }

    void collapse_halfstep(int k)
    {
        for (int i = 0; i < max_x; i++)
        {
            for (int j = 0; j < max_y; j++)
            {

                U_mut(i, j, k) = U(i, j, k + 1);
                U_mut(i, j, k + 1) = 0;
            }
        }
    }
};

double shmain(int w, int h, int d)
{
    griddy x = griddy(w, h, d);
    x.annul();
    x.set_start();
    for (int k = 1; k < d; k++)
    {
        x.full_halfstep_unclear_x(k);
        x.full_halfstep_unclear_y(k + 1);
        x.collapse_halfstep(k);
        x.set_start();
    }
    griddy y = griddy(x);
    y.cheat_set_correct();
    griddy z = x.diff(y);
    return z.square_error();
}

int main()
{
    printf("\033[2J");
    int w = 10;
    int h = 10;
    int d = 10;
    double er = shmain(w, h, d);
    printf("error is %lf\n", er);
    return 0;
}