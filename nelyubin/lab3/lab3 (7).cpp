#include <iostream>
#include <cmath>
#include <vector>

const double MAX_X = M_PI / 2;
const double MAX_Y = M_PI / 2;

double phi_X0(double y)
{
    return exp(-y) * cos(y);
}

double phi_XL(double y)
{
    return 0;
}

double phi_Y0(double x)
{
    return cos(x);
}

double phi_YL(double t)
{
    return 0;
}

double correct(double x, double y)
{
    return exp(-y) * cos(x) * cos(y);
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

/// DECISION: i w IS X AND j h IS Y
class griddy
{
    int w;
    int h;
    double hx;
    double hy;
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
        hx = MAX_X / (w - 1);
        hy = MAX_Y / (h - 1);
        vals = std::vector<double>(w * h);
    }
    griddy(griddy &other)
    {
        w = other.w;
        h = other.h;
        hx = other.hx;
        hy = other.hy;
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
        return pow(1 / hx, 2) * hy;
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
                *U_mut(i, j) = correct(i * hx, j * hy);
            }
        }
    }
    void set_start()
    {

        for (int i = 0; i < w; i++)
        {
            *U_mut(i, 0) = phi_Y0(i * hx);
            *U_mut(i, h - 1) = phi_YL(i * hx);
        }
        for (int j = 0; j < h; j++)
        {
            *U_mut(0, j) = phi_X0(j * hy);
            *U_mut(w - 1, j) = phi_XL(j * hy);
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

    void lerp()
    {
        for (int i = 1; i < w - 1; i++)
        {
            double low_x = U(0, 0) * (w - i) / (w + 1);
            low_x += U(w - 1, 0) * (i + 1) / (w + 1);
            double high_x = U(0, h - 1) * (w - i) / (w + 1);
            high_x += U(w - 1, h - 1) * (i + 1) / (w + 1);
            for (int j = 1; j < h - 1; j++)
            {
                double val = low_x * (h - j) / (h + 1);
                val += high_x * (j + 1) / (h + 1);
                *U_mut(i, j) = val;
            }
        }
    }

    void liebman()
    {
        std::vector<double> new_vals(w * h);
        for (int i = 1; i < w - 1; i++)
        {
            for (int j = 1; j < h - 1; j++)
            {
                double p1 = hy * hy * (U(i + 1, j) + U(i - 1, j));
                double p2 = hx * hx * (U(i, j + 1) + U(i, j - 1));
                double p3 = hx * hx * hy * (U(i, j + 1) - U(i, j - 1));
                double p4k = 3 * hx * hx * hy * hy;
                double k1 = 2 * hy * hy;
                double k2 = 2 * hx * hx;
                double all_k = k1 + k2;
                double d = p1 + p2 + p3 + p4k * U(i, j);
                new_vals[i * w + j] = d / all_k;
            }
        }
        for (int i = 1; i < w - 1; i++)
        {
            for (int j = 1; j < h - 1; j++)
            {
                *U_mut(i, j) = new_vals[i * w + j];
            }
        }
    }

    void zeidel()
    {
        for (int i = 1; i < w - 1; i++)
        {
            for (int j = 1; j < h - 1; j++)
            {
                double p1 = hy * hy * (U(i + 1, j) + U(i - 1, j));
                double p2 = hx * hx * (U(i, j + 1) + U(i, j - 1));
                double p3 = hx * hx * hy * (U(i, j + 1) - U(i, j - 1));
                double p4k = 3 * hx * hx * hy * hy;
                double k1 = 2 * hy * hy;
                double k2 = 2 * hx * hx;
                double all_k = k1 + k2;
                double d = p1 + p2 + p3 + p4k * U(i, j);
                *U_mut(i, j) = d / all_k;
            }
        }
    }
};

double const ERROR_MARGIN = 0.0001;

double shmain(int w, int h)
{
    griddy x = griddy(w, h);
    // printf("kurant %lf\n", x.krant());
    x.annul();
    // *x.U_mut(0, 0) = -1;
    // *x.U_mut(w - 1, h - 1) = 1;
    // x.cheat_set_correct();
    x.set_start();
    x.lerp();
    int j;
    for (j = 1; j < 100; j++)
    {
        griddy past_x = griddy(x);
        // x.liebman();
        x.zeidel();
        past_x = past_x.diff(x);
        if (past_x.square_error() < ERROR_MARGIN)
        {
            break;
        }
    }
    printf("took %d iterations\n", j);
    griddy y = griddy(x);
    y.cheat_set_correct();
    griddy z = x.diff(y);
    // z.print();
    // for (int j = 0; j < h; j++)
    // {
    //     printf("%lf\n", z.n_row_error(j));
    // }
    return z.square_error();
}

int main()
{
    int w = 10;
    int h = 10;
    double er = shmain(w, h);
    printf("error is %lf\n", er);
    return 0;
}