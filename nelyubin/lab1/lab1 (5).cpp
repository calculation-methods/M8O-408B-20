#include <iostream>
#include <cmath>
#include <vector>

// dU = a d2U + b dU
// dt     dx2 +   dx
//  a>0 b>0
//  Ux(0,t)-U(0,t)=-exp(-at)(cos(bt)+sin(bt))
//  Ux(pi,t)-U(pi,t)=exp(-at)(cos(bt)+sin(bt))
//  U(x,0)=cos(x)

// dU(0,t)/dx = U(0,t)-exp(at)(cos(bt)+sin(bt))

const double MAX_X = M_PI;
const double MAX_T = 1;
/// in methods its a^2
const double A = 0.1;
const double B = 0.1;
/// unused
const double C = 0;

const double ALPHA = 1;
const double BETA = -1;
const double GAMMA = 1;
const double DELTA = -1;

double initial(double x)
{
    return cos(x);
}

double phi0(double t)
{
    return -exp(-A * t) * (cos(B * t) + sin(B * t));
}

double phiL(double t)
{
    return exp(-A * t) * (cos(B * t) + sin(B * t));
}

// U(x,t)=exp(-at)cos(x+bt)

double correct(double x, double t)
{
    return exp(-A * t) * cos(x + B * t);
}

// Ux = dU/dx

// todo
//  еонучно ращностная схема, явная-неявная; Схема Кранка Николсана
//  граничные условия апроксимация:
//       двухточечная с 1 порядком
//       трёхточечная 2 поядок
//       двухточечная со вторым порядком

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

enum border_precision
{
    point_2_order_1,
    point_2_order_2,
    point_3_order_2
};

/// DECISION: i w IS SPACE AND j h IS TIME
class griddy
{
    int w;
    int h;
    border_precision prec;
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
    void set_prec(border_precision x)
    {
        prec = x;
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
        return pow(1 / hsh, 2) * tau * A;
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
            for (int j = 1; j < h; j++)
            {
                *U_mut(i, j) = correct(i * hsh, j * tau);
            }
        }
    }
    void set_start()
    {
        annul();
        for (int i = 0; i < w; i++)
        {
            *U_mut(i, 0) = initial(i * hsh);
            for (int j = 1; j < h; j++)
            {
                *U_mut(i, j) = 0;
            }
        }
    }

    void clear_2p_1o(int k)
    {
        *U_mut(0, k + 1) = -(ALPHA / hsh) * U(1, k + 1) + phi0(k + 1);
        *U_mut(0, k + 1) /= (BETA - ALPHA / hsh);
        *U_mut(w - 1, k + 1) = (GAMMA / hsh) / (DELTA + GAMMA / hsh) * U(w - 2, k + 1) + phiL(k + 1) / (DELTA + GAMMA / hsh);
    }

    void clear_3p_2o(int k)
    {
        *U_mut(0, k + 1) = phi0(k + 1) - ALPHA * (4 * U(1, k + 1) - U(2, k + 1)) / (2 * hsh);
        *U_mut(0, k + 1) /= BETA - 3 * ALPHA / (2 * hsh);
        *U_mut(w - 1, k + 1) = phiL(k + 1) - GAMMA * (U(w - 3, k + 1) - 4 * U(w - 2, k + 1)) / (2 * hsh);
        *U_mut(w - 1, k + 1) /= DELTA + 3 * GAMMA / (2 * hsh);
    }

    void clear_2p_2o(int k)
    {
        double a0 = 0;
        double b0 = (2 * A) / hsh + hsh / tau - C * hsh - (BETA / ALPHA) * (2 * A - B * hsh);
        double c0 = -2 * A / hsh;
        double d0 = hsh / tau * U(0, k) - phi0(k + 1) * (2 * A - B * hsh) / ALPHA;
        *U_mut(0, k + 1) = (d0 - c0 * U(1, k + 1)) / b0;

        double aN = -2 * A / hsh;
        double bN = (2 * A) / hsh + hsh / tau - C * hsh - (DELTA / GAMMA) * (2 * A - B * hsh);
        double cN = 0;
        double dN = hsh / tau * U(w - 1, k) - phiL(k + 1) * (2 * A - B * hsh) / GAMMA;
        *U_mut(w - 1, k + 1) = (dN - cN * U(w - 2, k + 1)) / bN;
    }

    void clear_scheme(int k)
    {
        for (int i = 1; i < w - 1; i++)
        {
            *U_mut(i, k + 1) = krant() * (U(i - 1, k) + U(i + 1, k)) + (1 - 2 * krant()) * U(i, k);
            *U_mut(i, k + 1) += tau * B / (2 * hsh) * (U(i + 1, k) - U(i - 1, k));
        }
        switch (prec)
        {
        case point_2_order_1:
            clear_2p_1o(k);
            break;
        case point_2_order_2:
            clear_2p_2o(k);
            break;
        case point_3_order_2:
            clear_3p_2o(k);
            break;
        }
    }

    void unclear_2p_1o(std::vector<double> &ar, std::vector<double> &br, std::vector<double> &cr, std::vector<double> &dr, int k)
    {
        for (int i = 0; i < w; i++)
        {
            ar[i] = krant();
            br[i] = -(1 + 2 * krant());
            cr[i] = krant();
            dr[i] = -U(i, k);
        }
        ar[0] = 0;
        br[0] = BETA - ALPHA / hsh;
        cr[0] = ALPHA / hsh;
        dr[0] = phi0(k + 1) / (BETA - ALPHA / hsh);
        int l = w - 1;
        ar[l] = -GAMMA / hsh;
        br[l] = DELTA + GAMMA / hsh;
        cr[l] = 0;
        dr[l] = phiL(k + 1) / (DELTA + GAMMA / hsh);
    }

    void unclear_3p_2o(std::vector<double> &ar, std::vector<double> &br, std::vector<double> &cr, std::vector<double> &dr, int k)
    {
        for (int i = 0; i < w; i++)
        {
            ar[i] = krant();
            br[i] = -(1 + 2 * krant());
            cr[i] = krant();
            dr[i] = -U(i, k);
        }

        ar[0] = 0;
        br[0] = BETA - 3 * ALPHA / (2 * hsh);
        cr[0] = -2 * ALPHA / hsh;
        double e0 = ALPHA / (2 * hsh);
        dr[0] = phi0(k + 1) / br[0];
        double k0 = e0 / cr[1];
        br[0] -= ar[1] * k0;
        cr[0] -= br[1] * k0;
        e0 -= cr[1] * k0;
        if (e0 != 0)
        {
            printf("ERROROROR 0 %lf\n", e0);
        }
        dr[0] -= dr[1] * k0;

        ar[w - 1] = 2 * GAMMA / hsh;
        br[w - 1] = DELTA + 3 * GAMMA / (2 * hsh);
        cr[w - 1] = 0;
        double el = -GAMMA / (2 * hsh);
        dr[w - 1] = phiL(k + 1) / br[w - 1];
        double kl = el / cr[w - 2];
        br[w - 1] -= cr[w - 2] * kl;
        ar[w - 1] -= br[w - 2] * kl;
        el -= ar[w - 2] * kl;
        if (el != 0)
        {
            printf("ERROROROR L %lf\n", e0);
        }
        dr[w - 1] -= dr[w - 2] * kl;
    }

    void unclear_2p_2o(std::vector<double> &ar, std::vector<double> &br, std::vector<double> &cr, std::vector<double> &dr, int k)
    {
        for (int i = 0; i < w; i++)
        {
            ar[i] = -(A / pow(hsh, 2) - B / (2 * hsh));
            br[i] = 2 * A / pow(hsh, 2) + 1 / tau - C;
            cr[i] = -(A / pow(hsh, 2) + B / (2 * hsh));
            dr[i] = U(i, k) / tau;
        }
        ar[0] = 0;
        cr[w - 1] = 0;
    }

    void unclear_scheme(int k)
    {
        // this is 2-2
        std::vector<double> ar(w, 0);
        std::vector<double> br(w, 0);
        std::vector<double> cr(w, 0);
        std::vector<double> dr(w, 0);
        std::vector<double> xr(w, 0);

        switch (prec)
        {
        case point_2_order_1:
            unclear_2p_1o(ar, br, cr, dr, k);
            break;
        case point_2_order_2:
            unclear_2p_2o(ar, br, cr, dr, k);
            break;
        case point_3_order_2:
            unclear_3p_2o(ar, br, cr, dr, k);
            break;
        }

        // for (int i = 0; i < w; i++)
        // {
        //     printf("a=%2.3lf  a=%2.3lf  a=%2.3lf  a=%2.3lf  \n", ar[i], br[i], cr[i], dr[i]);
        // }

        xr = progon(ar, br, cr, dr);
        for (int i = 0; i < w; i++)
        {
            *U_mut(i, k + 1) = xr[i];
        }
    }

    void krank_nickolson(int k)
    {
        double theta = 0.5;
        std::vector<double> half(w, 0);
        clear_scheme(k);
        for (int i = 0; i < w; i++)
        {
            half[i] = U(i, k);
        }
        unclear_scheme(k);
        for (int i = 0; i < w; i++)
        {
            *U_mut(i, k) = U(i, k) * (1 - theta) + half[i] * theta;
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
    x.set_prec(point_3_order_2);
    x.set_start();
    for (int j = 0; j < h - 1; j++)
    {
        // x.clear_scheme(j);
        x.unclear_scheme(j);
        // x.krank_nickolson(j);
    }
    griddy y = griddy(x);
    y.cheat_set_correct();
    griddy z = x.diff(y);
    z.print();
    // for (int j = 0; j < h; j++)
    // {
    //     printf("%lf\n", z.n_row_error(j));
    // }
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