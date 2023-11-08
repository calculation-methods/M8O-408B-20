#include <iostream>
#include <cmath>
#include <vector>

const double MAX_X = M_PI;
const double MAX_T = 1;

double f(double x, double t)
{
    return sin(x) * exp(-t);
}

double psi1(double x)
{
    return cos(x);
}

double psi2(double x)
{
    return -cos(x);
}

double phi0(double t)
{
    return exp(-t);
}

double phiL(double t)
{
    return -exp(-t);
}

double correct(double x, double t)
{
    return exp(-t) * cos(x);
}

std::vector<double> progon(std::vector<double> ar, std::vector<double> br, std::vector<double> cr, std::vector<double> dr)
{
    // printf("%d %d %d %d\n", ar.size(), br.size(), cr.size(), dr.size());
    int len = br.size();
    if (dr.size() != len)
    {
        throw std::invalid_argument("bad array dimensions");
    }
    if (ar.size() == len)
    {
        if (ar[0] != 0)
        {
            throw std::invalid_argument("a[0] must be 0");
        }
    }
    else if (ar.size() + 1 == len)
    {
        ar.insert(ar.begin(), 0);
    }
    else
    {
        throw std::invalid_argument("bad array dimensions");
    }
    if (cr.size() == len)
    {
        if (cr[len - 1] != 0)
        {
            throw std::invalid_argument("c[LAST] must be 0");
        }
    }
    else if (cr.size() + 1 == len)
    {
        cr.insert(ar.end(), 0);
    }
    else
    {
        throw std::invalid_argument("bad array dimensions");
    }
    std::vector<double> P;
    std::vector<double> Q;
    for (int i = 0; i < len; i++)
    {
        double help = br[i];
        if (i > 0)
        {
            help += ar[i] * P[i - 1];
        }
        P.push_back(-cr[i] / help);
        double help2 = dr[i];
        if (i > 0)
        {
            help2 -= ar[i] * Q[i - 1];
        }
        Q.push_back(help2 / help);
    }

    // for (int i = 0; i < len; i++)
    // {
    //     printf("P[%d] = %lf\tQ[%d] = %lf\n", i, P[i], i, Q[i]);
    // }

    std::vector<double> result(len, Q[len - 1]);
    for (int i = len - 2; i >= 0; i--)
    {
        result[i] = P[i] * result[i + 1] + Q[i];
    }
    return result;
}

/// DECISION: i w IS SPACE AND j h IS TIME
class griddy
{
    int w;
    int h;
    double hsh;
    double tau;
    std::vector<double> vals;

    void process_arguments(int i, int j)
    {
        if ((i < 0) || (i >= w))
        {
            throw std::invalid_argument("index i outside of range");
        }
        if ((j < 0) || (j >= h))
        {
            throw std::invalid_argument("index j outside of range");
        }
    }

public:
    griddy(int width, int height)
    {
        w = width;
        h = height;
        hsh = MAX_X / (w - 1);
        tau = MAX_T / (h - 1);
        vals = std::vector<double>(w * h);
    }
    griddy(griddy &other)
    {
        w = other.w;
        h = other.h;
        hsh = other.hsh;
        tau = other.tau;
        vals = std::vector<double>(w * h);
        for (int i = 0; i < w * h; i++)
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
    double krant()
    {
        return pow(1 / hsh, 2) * tau;
    }
    double f_proc(int i, int j)
    {
        return f(i * hsh, j * tau);
    }
    double *U_mut(int i, int j)
    {
        process_arguments(i, j);
        int idx = i * h + j;
        return &vals[idx];
    }
    double U(int i, int j)
    {
        return *U_mut(i, j);
    }
    void print()
    {
        // for (int i = 0; i < vals.size(); i++)
        // {
        //     printf("%3.0lf ", vals[i]);
        // }
        // printf("\n");

        // for (int i = 0; i < w; i++)
        // {
        for (int j = 0; j < h; j++)
        {
            for (int i = 0; i < w; i++)
            {
                printf("%8.4lf ", U(i, j));
            }
            printf("\n");
        }
    }
    void cheat_set_correct()
    {
        for (int i = 0; i < w; i++)
        {
            for (int j = 0; j < h; j++)
            {
                *U_mut(i, j) = correct(i * hsh, j * tau);
            }
        }
    }
    void set_start()
    {
        for (int i = 0; i < w; i++)
        {
            *U_mut(i, 0) = psi1(i * hsh);
        }
        for (int j = 1; j < h; j++)
        {
            *U_mut(0, j) = phi0(j * tau);
            *U_mut(w - 1, j) = phiL(j * tau);
        }
    }

    void second_level_1o()
    {
        for (int i = 0; i < w; i++)
        {
            *U_mut(i, 1) = U(i, 0) + tau * psi2(i * hsh);
        }
    }

    void second_level_2o()
    {
        for (int i = 0; i < w; i++)
        {
            *U_mut(i, 1) = psi1(i * hsh) + tau * psi2(i * hsh) + psi1(i * hsh) * tau * tau / 2;
        }
    }

    void cross(int j)
    {
        for (int i = 1; i < w - 2; i++)
        {
            double p = (1 + hsh) * U(i + 1, j) - (2 + hsh) * U(i, j) + U(i - 1, j);
            p /= hsh * hsh;
            p += U(i, j) + f_proc(i, j);
            p *= tau * tau;
            p += (2 + 3 * tau) * U(i, j) - U(i, j - 1);
            p /= 3 * tau + 1;
            *U_mut(i, j + 1) = p;
        }
    }

    void kris_cross(int j)
    {
        for (int i = 1; i < w - 2; i++)
        {
            double p = (2 + hsh) * U(i + 1, j) - 4 * U(i, j) + (2 - hsh) * U(i - 1, j);
            p /= 2 * hsh * hsh;
            p += U(i, j) + f_proc(i, j);
            p *= 2 * tau * tau;
            p -= -4 * U(i, j) + (2 - 3 * tau) * U(i, j - 1);
            p /= 2 + 3 * tau;
            *U_mut(i, j + 1) = p;
        }
    }

    void big_ol_t(int j)
    {
        std::vector<double> as(0);
        std::vector<double> bs(0);
        std::vector<double> cs(0);
        std::vector<double> ds(0);
        for (int i = 1; i < w - 2; i++)
        {
            double sa = -(2 - tau) / (2 * hsh * hsh);
            double sb = -(-4) / (2 * hsh * hsh);
            double sc = -(2 + tau) / (2 * hsh * hsh);
            sb += (2 + 3 * tau) / (2 * tau * tau);
            double sd = U(i, j) + f_proc(i, j);
            sd -= (-4 * U(i, j) + (2 - 3 * tau) * U(i, j - 1)) / (2 * tau * tau);
            as.push_back(sa);
            bs.push_back(sb);
            cs.push_back(sc);
            ds.push_back(sd);
        }
        ds[0] -= as[0] * U(0, j + 1);
        as[0] = 0;
        int l = cs.size() - 1;
        ds[l] -= cs[l] * U(w - 1, j + 1);
        cs[l] = 0;
        // for (int i = 0; i < as.size(); i++)
        // {
        //     for (int j = 0; j < i; j++)
        //     {
        //         printf("            ");
        //     }
        //     printf("%7.3lf a + %7.3lf b + %7.3lf c", as[i], bs[i], cs[i]);
        //     for (int j = i + 1; j < as.size(); j++)
        //     {
        //         printf("            ");
        //     }
        //     printf(" = %7.3lf d\n", ds[i]);
        // }
        std::vector<double> res = progon(as, bs, cs, ds);
        for (int i = 1; i < w - 2; i++)
        {
            *U_mut(i, j + 1) = res[i];
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
    double n_row_error(int j)
    {
        double res = 0;
        for (int i = 0; i < w; i++)
        {
            res += pow(U(i, j), 2);
        }
        return sqrt(res);
    }
    griddy diff(griddy &other)
    {
        if ((other.w != w) || (other.h != h))
        {
            throw std::invalid_argument("dimensions do not allign");
        }
        griddy rez = griddy(w, h);
        for (int i = 0; i < w; i++)
        {
            for (int j = 0; j < h; j++)
            {
                *rez.U_mut(i, j) = U(i, j) - other.U(i, j);
            }
        }
        return rez;
    }
};

double shmain(int w, int h)
{
    griddy x = griddy(w, h);
    printf("kurant %lf\n", x.krant());
    x.annul();
    x.cheat_set_correct();
    x.set_start();
    x.second_level_2o();
    for (int j = 1; j < h - 1; j++)
    {
        // x.kris_cross(j);
        x.big_ol_t(j);
        // break;
    }
    griddy y = griddy(x);
    y.cheat_set_correct();
    griddy z = x.diff(y);
    // z.print();
    for (int j = 0; j < h; j++)
    {
        printf("%lf\n", z.n_row_error(j));
    }
    return z.n_row_error(h - 1);
}

int main()
{
    int w = 10;
    int h = 7;
    double er = shmain(w, h);
    printf("%lf\n", er);
    return 0;
}