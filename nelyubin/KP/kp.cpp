#include <iostream>
#include <vector>

typedef double db;
typedef std::vector<db> da;
typedef size_t st;

class mtrx
{
public:
    st w, h;
    std::vector<st> idx;
    da vl;
    mtrx(st ww, st hh, bool warn = false)
    {
        if (warn)
        {
            if (hh > ww)
            {
                std::cout << "excess lines\n";
            }
            if (hh < ww)
            {
                std::cout << "infinite solutions\n";
            }
        }
        w = ww;
        h = hh;
        idx.clear();
        vl.clear();
    }
    mtrx(const mtrx &tar)
    {
        w = tar.w;
        h = tar.h;
        idx = tar.idx;
        vl = tar.vl;
    }

    mtrx(const da &tar, bool T)
    {
        if (T)
        {
            w = 1;
            h = tar.size();
        }
        else
        {
            w = tar.size();
            h = 1;
        }
        idx.clear();
        for (st i = 0; i < tar.size(); i++)
        {
            idx.push_back(i);
        }
        vl = tar;
    }

    st d2s(st i, st j) const
    {
        st z = i * w + j;
        int l = 0;
        int r = vl.size();
        int m;
        while (l + 1 < r)
        {
            m = (l + r) / 2;
            if (idx[m] == z)
            {
                return (st)m;
            }
            if (idx[m] > z)
            {
                r = m;
                continue;
            }
            if (idx[m] < z)
            {
                l = m;
                continue;
            }
        }
        return (st)l;
    }

    db &operator()(st ii, st jj)
    {
        if (idx.size() == 0)
        {
            idx.push_back(ii * w + jj);
            vl.push_back(0);
            return vl[0];
        }
        st ix = d2s(ii, jj);
        if ((ix >= idx.size()) || (idx[ix] != ii * w + jj))
        {
            idx.insert(idx.begin() + ix + 1, 1, ii * w + jj);
            vl.insert(vl.begin() + ix + 1, 1, 0);
            return vl[ix + 1];
        }
        return vl[ix];
    }

    db get(st ii, st jj) const
    {
        if (idx.size() == 0)
        {
            return 0;
        }
        st ix = d2s(ii, jj);
        if ((ix >= idx.size()) || (idx[ix] != ii * w + jj))
        {
            return 0;
        }
        return vl[ix];
    }

    std::pair<da, da> randomFill(st chance, st range, bool tr)
    {
        st rnd;
        std::pair<da, da> rez;
        for (st i = 0; i < h; i++)
        {
            (*this)(i, i) = ((double)rand() / RAND_MAX * 2 - 1) * range;
            for (st j = (i + 1) * tr; j < w; j++)
            {
                rnd = (double)rand() / RAND_MAX * 100;
                if (rnd < chance)
                {
                    (*this)(i, j) = (double)(rand() % range) * 2 - range;
                }
            }
        }
        for (st i = 0; i < h; i++)
        { // these are X
            rez.first.push_back(rand() - 5);
        }
        for (st i = 0; i < h; i++)
        { // these are B
            rez.second.push_back(0);
            for (st j = 0; j < h; j++)
            {
                rez.second[i] += this->get(i, j) * rez.first[j];
            }
        }
        return rez;
    }

    void print() const
    {
        st k = 0;
        for (st i = 0; i < h; i++)
        {
            for (st j = 0; j < w; j++)
            {
                if ((k < idx.size()) && (idx[k] == i * w + j))
                {
                    std::cout << vl[k] << "\t";
                    k++;
                }
                else
                {
                    std::cout << "-\t";
                }
            }
            std::cout << "\n";
        }
    }

    void print(const mtrx &tar) const
    {
        st k = 0;
        st ks = 0;
        for (st i = 0; i < h; i++)
        {
            for (st j = 0; j < w; j++)
            {
                if ((k < idx.size()) && (idx[k] == i * w + j))
                {
                    std::cout << vl[k] << "\t";
                    k++;
                }
                else
                {
                    std::cout << "-\t";
                }
            }
            std::cout << "\t\t";
            for (st j = 0; j < tar.w; j++)
            {
                if ((ks < tar.idx.size()) && (tar.idx[ks] == i * tar.w + j))
                {
                    std::cout << tar.vl[ks] << "\t";
                    ks++;
                }
                else
                {
                    std::cout << "-\t";
                }
            }
            std::cout << "\n";
        }
    }

    void clear()
    {
        for (st i = 0; i < idx.size(); i++)
        {
            if (vl[i] == 0)
            {
                idx.erase(idx.begin() + i);
                vl.erase(vl.begin() + i);
                i--;
            }
        }
    }

    void e()
    {
        for (st i = 0; i < w; i++)
        {
            for (st j = 0; j < h; j++)
            {
                (*this)(i, j) = (int)(i == j);
            }
        }
        this->clear();
    }

    mtrx t()
    {
        db hl;
        mtrx rez(h, w);
        for (st i = 0; i < h; i++)
        {
            for (st j = 0; j < w; j++)
            {
                hl = this->get(i, j);
                if (hl != 0)
                {
                    rez(j, i) = hl;
                }
            }
        }
        return rez;
    }

    mtrx operator*(mtrx &b)
    {
        int a1, ab, b2;
        a1 = h;
        ab = b.h;
        b2 = b.w;
        if (w != ab)
        {
            std::cout << "cannot multiply matrixes\n";
        }
        mtrx rez(b2, a1);
        da col;
        col.clear();
        double hh = 0;
        for (int i = 0; i < a1; i++)
        {
            for (int j = 0; j < b2; j++)
            {
                for (int k = 0; k < ab; k++)
                {
                    hh += this->get(i, k) * b.get(k, j);
                }
                rez(i, j) = hh;
                hh = 0;
            }
        }
        rez.clear();
        return rez;
    }

    mtrx operator*(double a)
    {
        mtrx rez(*this);
        for (st i = 0; i < rez.vl.size(); i++)
        {
            rez.vl[i] *= a;
        }
        rez.clear();
        return rez;
    }

    mtrx operator+(const mtrx &b)
    {
        if ((w != b.w) || (h != b.h))
        {
            std::cout << "cannot sum matrixes\n";
        }
        db hl = 0;
        mtrx rez(*this);
        for (st i = 0; i < h; i++)
        {
            for (st j = 0; j < w; j++)
            {
                hl = b.get(i, j);
                if (hl != 0)
                {
                    rez(i, j) += hl;
                }
            }
        }
        rez.clear();
        return rez;
    }

    mtrx operator-(const mtrx &b)
    {
        if ((w != b.w) || (h != b.h))
        {
            std::cout << "cannot sub matrixes\n";
        }
        db hl = 0;
        mtrx rez(*this);
        for (st i = 0; i < h; i++)
        {
            for (st j = 0; j < w; j++)
            {
                hl = b.get(i, j);
                if (hl != 0)
                {
                    rez(i, j) -= hl;
                }
            }
        }
        rez.clear();
        return rez;
    }

    bool operator==(const mtrx &b)
    {
        if (w != b.w)
        {
            return false;
        }
        if (h != b.h)
        {
            return false;
        }
        if (idx.size() != b.idx.size())
        {
            return false;
        }
        for (int i = 0; i < idx.size(); i++)
        {
            if (idx[i] != b.idx[i])
            {
                return false;
            }
            if (vl[i] != b.vl[i])
            {
                return false;
            }
        }
        return true;
    }
};

mtrx operator*(const double &a, const mtrx &tar)
{
    mtrx rez(tar);
    for (st i = 0; i < rez.vl.size(); i++)
    {
        rez.vl[i] *= a;
    }
    rez.clear();
    return rez;
}

db dot(const mtrx &a, const mtrx &b)
{
    db rez = 0;
    if (a.h == 1)
    {
        if (b.h == 1)
        {
            if (a.w != b.w)
            {
                std::cout << "cant dot 1\n";
            }
            for (st i = 0; i < a.w; i++)
            {
                rez += a.get(i, 0) * b.get(i, 0);
            }
        }
        else if (b.w == 1)
        {
            if (a.w != b.h)
            {
                std::cout << "cant dot 2\n";
            }
            for (st i = 0; i < a.w; i++)
            {
                rez += a.get(i, 0) * b.get(0, i);
            }
        }
    }
    else if (a.w == 1)
    {
        if (b.h == 1)
        {
            if (a.h != b.w)
            {
                std::cout << "cant dot 3\n";
            }
            for (st i = 0; i < a.h; i++)
            {
                rez += a.get(0, i) * b.get(i, 0);
            }
        }
        else if (b.w == 1)
        {
            if (a.h != b.h)
            {
                std::cout << "cant dot 3\n";
            }
            for (st i = 0; i < a.h; i++)
            {
                rez += a.get(0, i) * b.get(0, i);
            }
        }
    }
    return rez;
}

mtrx multMatrix(mtrx &a, mtrx &b)
{
    int a1, ab, b2;
    a1 = a.h;
    ab = b.h;
    b2 = b.w;
    if (a.w != ab)
    {
        std::cout << "cannot multiply matrixes\n";
    }
    mtrx rez(b2, a1);
    da col;
    col.clear();
    double h = 0;
    for (int i = 0; i < a1; i++)
    {
        for (int j = 0; j < b2; j++)
        {
            for (int k = 0; k < ab; k++)
            {
                h += a.get(i, k) * b.get(k, j);
            }
            rez(i, j) = h;
            h = 0;
        }
    }
    rez.clear();
    return rez;
}

int main()
{
    db er = 1.0 / 1000000;
    std::cout.precision(3);
    srand(time(0));
    std::pair<da, da> data;
    mtrx m(10, 10), ep(m);
    data = m.randomFill(5, 10, true);
    ep.randomFill(5, 5, false);
    mtrx B(data.second, true);
    mtrx ans(data.first, true);
    mtrx A = multMatrix(ep, m);
    mtrx sb = multMatrix(ep, B);
    std::cout << "\n";
    A.print(sb);
    std::cout << "\n";

    mtrx x0 = sb * 0;

    mtrx r0 = sb - (A * x0);
    mtrx v0 = r0 * 0.0;
    mtrx p0 = v0;
    db ro0 = 1, al = 1, w0 = 1;
    db nroi, bt, oroi = ro0;
    db owi = w0, nwi = w0;
    mtrx ori = r0, nri = r0;
    mtrx npi = ori, opi = ori;
    mtrx ovi = v0, nvi = v0;
    mtrx h = r0, s = r0, t = r0;
    mtrx oxi = x0, nxi = x0;
    for (st i = 0; i < 160; i++)
    {
        nroi = dot(r0, ori);
        bt = (nroi / oroi) * (al / owi);
        npi = ori + bt * (opi - owi * ovi);
        nvi = A * npi;
        al = nroi / dot(r0, nvi);
        h = oxi + al * npi;
        s = ori - al * nvi;
        t = A * s;
        nwi = dot(t, s) / dot(t, t);
        nxi = h + nwi * s;
        nri = s - nwi * t;
        if (nxi == oxi)
        {
            break;
        }
        oxi = nxi;
        oroi = nroi;
        opi = npi;
        ori = nri;
        owi = nwi;
        ovi = nvi;
    }
    nxi.print(ans);

    return 0;
}