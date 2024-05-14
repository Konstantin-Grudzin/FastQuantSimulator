#define _CRT_SECURE_NO_WARNINGS
#pragma once
//#pragma omp parallel for
//#pragma omp parallel
#include <iostream>
#include<string>
#include <vector>
#include<random>
#include<ctime>
#include<bitset>
#include<string>
#include<unordered_map>
#include<cmath>
#include<iomanip>

using namespace std;

typedef long double ld;
static ld aa = 1 / sqrt(2);
static ld pi = acos(-1);
static ld eps = 1e-15;

static int gcd(int a, int b)
{
    if (b == 0)
        return a;
    return gcd(b, a % b);
}

static int RevModNum(int a, int N)
{
    if (a == 1)
        return 1;
    int q0, q1, q2, r, s1, s2, t1, t2, tmp, res;
    t2 = s1 = 1;
    t1 = s2 = 0;
    r = q1 = N;
    q2 = a;

    while (r > 0)
    {
        q0 = q1;
        tmp = q2;
        q2 = q1 / q2;
        q1 = tmp;
        r = q0 - q1 * q2;

        tmp = s1;
        s1 = s2;
        s2 = tmp - s1 * q2;

        tmp = t1;
        t1 = t2;
        t2 = tmp - t1 * q2;

        q2 = r;

        if (r)
            res = t2;
    }
    while (res < 0)
        res += N;
    return res;
}

static string toBits(int i, int qc)
{
    string ans = "";
    for (int j = 0; j < qc; ++j)
    {
        ans = to_string(i % 2) + ans;
        i /= 2;
    }
    return ans;
}

static int pow2(int a, int n, int mod)
{
    int res = 1;
    while (n) {
        if (n & 1)
            res = (a * res % mod + mod) % mod;
        a = (a * a % mod + mod) % mod;
        n >>= 1;
    }
    return res;
}

static int pow2(int a, int n)
{
    int res = 1;
    while (n) {
        if (n & 1)
            res = a * res;
        a = a * a;
        n >>= 1;
    }
    return res;
}

class complex
{
    ld a;
    ld b;
public:
    complex(ld _a = 0, ld _b = 0) : a(_a), b(_b)
    {

    }
    complex operator+(complex y)
    {
        complex z(a + y.a, b + y.b);
        return z;
    }
    complex operator-(complex y)
    {
        complex z(a - y.a, b - y.b);
        return z;
    }
    complex operator*(complex y)
    {
        complex z(a * y.a - b * y.b, a * y.b + b * y.a);
        return z;
    }
    complex operator/ (ld x)
    {
        complex z(a / x, b / x);
        return z;
    }
    complex operator* (ld x)
    {
        complex z(a * x, b * x);
        return z;
    }

    bool operator==(complex v)
    {
        return (a == v.a) && (b == v.b);
    }
    bool operator!=(complex v)
    {
        return !(*this == v);
    }
    ld SqMod()
    {
        return a * a + b * b;
    }

    ld SqrtMod()
    {
        return sqrt(a * a + b * b);
    }

    friend istream& operator>>(istream& in, complex a)
    {
        in >> a.a >> a.b;
        return in;
    }

    friend ostream& operator<<(ostream& out, complex a)
    {
        out << a.a;
        if (a.b <= eps)
        {
            return out;
        }

        if (a.b > 0)
            out << " +" << a.b << " * i";
        else
            out << " " << a.b << " * i";
        return out;
    }

};

static complex exp_i_phi(double phi)
{
    complex a(cos(phi), sin(phi));
    return a;
}

static complex obrk2(aa);

//всегда начинаем из 0000
class qMatrix
{
    int size_in_bits;
    atomic<int> GateCount = 0;
    vector<complex> v;
    ld error;
public:
    qMatrix(int size = 1,ld err = 0.001) : size_in_bits(size),error(err)
    {
        v.resize(1 << (size_in_bits));
        v[0] = 1;
    }

    inline int GetGateCount()
    {
        return GateCount;
    }

    complex& operator[](int ind)
    {
        return v[ind];
    }

    void reset()
    {
        GateCount = 0;
        for (int i = 0; i < v.size(); ++i)
        {
            if (i == 0) v[i] = 1;
            else v[i] = 0;
        }
    }

    const qMatrix& operator=(const qMatrix& A)
    {
        if (size_in_bits != A.size_in_bits)
        {
            size_in_bits = A.size_in_bits;

        }
        v.assign(begin(A.v), end(A.v));

        return *this;
    }

    void X(int ind)
    {
        GateCount++;
        int size = (1 << size_in_bits);
        for (int i = 0; i < size; ++i)
        {
            if (i & (1 << ind))
            {
                complex t = v[i ^ (1 << ind)];
                v[i ^ (1 << ind)] = v[i];
                v[i] = t;
            }
        }
    }


    void H(int ind)
    {
        GateCount++;
        qMatrix buff(size_in_bits);
        buff[0] = 0;
        int size = (1 << size_in_bits);
        for (int i = 0; i < size; ++i)
        {
            if (i & (1 << ind))
            {
                buff[i ^ (1 << ind)] = buff[i ^ (1 << ind)] + v[i] * obrk2;
                buff[i] = buff[i] - v[i] * obrk2;
            }
            else
            {
                buff[i ^ (1 << ind)] = buff[i ^ (1 << ind)] + v[i] * obrk2;
                buff[i] = buff[i] + v[i] * obrk2;
            }
        }
        *this = buff;
    }

    void CNOT(int cind, int chnind)
    {
        GateCount++;
        qMatrix a(size_in_bits);
        int size = 1 << (size_in_bits);
        for (int i = 0; i < size; ++i)
        {
            if (i & (1 << cind))
            {
                a[i ^ (1 << chnind)] = v[i];
            }
            else
            {
                a[i] = v[i];
            }
        }
        *this = a;
    }

    void CCNOT(int cind0, int cind1, int chnind)
    {
        GateCount++;
        qMatrix a(size_in_bits);
        int size = 1 << (size_in_bits);
        for (int i = 0; i < size; ++i)
        {
            if ((i & (1 << cind0)) && (i & (1 << cind1)))
            {
                a[i ^ (1 << chnind)] = v[i];
            }
            else
            {
                a[i] = v[i];
            }
        }
        *this = a;
    }

    void Ph(double theta, int bit)
    {
        GateCount++;
        std::random_device r;
        ld lower_bound = theta*(1-error);
        ld upper_bound = theta * (1 + error);
        std::uniform_real_distribution<ld> unif(lower_bound, upper_bound);
        std::default_random_engine re(r());
        ld alpha = unif(re);
        int size = 1 << (size_in_bits);
        complex exp = exp_i_phi(alpha);
#pragma omp parallel for
        for (int i = 0; i < size; ++i)
        {
            if ((i & (1 << bit)))
            {
                v[i] = v[i] * exp;
            }
        }
    }

    void R(int k, int bit)
    {
        ld x = pow2(2, k - 1);
        ld y = pi / x;
        Ph(y, bit);
    }

    void CPh(double theta, int cbit, int bit)
    {
        GateCount++;
        int size = 1 << (size_in_bits);
        complex exp = exp_i_phi(theta);
#pragma omp parallel for
        for (int i = 0; i < size; ++i)
        {
            if ((i & (1 << cbit)) && (i & (1 << bit)))
            {

                v[i] = v[i] * exp;
            }
        }
    }

    void CCPh(double theta, int cbit1, int cbit2, int bit)
    {
        GateCount++;
        int size = 1 << (size_in_bits);
        complex exp = exp_i_phi(theta);
#pragma omp parallel for
        for (int i = 0; i < size; ++i)
        {
            if ((i & (1 << cbit1)) && (i & (1 << cbit2)) && (i & (1 << bit)))
            {

                v[i] = v[i] * exp;
            }
        }
    }

    void QFT(int l, int r)
    {
        for (int i = r; i >= l; i--)
        {
            H(i);
            for (int j = i - 1, k = 4; j >= l; j--, k *= 2)
            {
                ld x = 2 * pi; x /= k;
                CPh(x, j, i);
            }
        }
    }

    void RQFT(int l, int r)
    {
        for (int i = l; i <= r; i++)
        {

            for (int j = i - 1, k = 4; j >= l; j--, k *= 2)
            {
                ld x = 2 * pi;
                x /= k;
                CPh(-x, j, i);
            }
            H(i);

        }
    }

    void FSUM(int l1, int r1, int l2, int r2)
    {

        for (int i = r2, t = 0; i >= l2; i--, t++)
        {
            for (int j = r1 - t, k = 2; j >= l1; j--, k *= 2)
            {
                ld x = 2 * pi; x /= k;
                CPh(x, j, i);
            }
        }
    }

    void FSUM(int c, int l1, int r1)
    {
        int l = 0;
        int r = r1 - l1;
        for (int i = r1, t = 0; i >= l1; i--, t++)
        {
            for (int j = r - t, k = 2; j >= l; j--, k *= 2)
            {
                if ((c >> j) & 1)
                {
                    ld x = 2 * pi; x /= k;
                    Ph(x, i);
                }
            }
        }
    }

    void CFSUM(int c, int cbit, int l1, int r1)
    {
        int l = 0;
        int r = r1 - l1;
        for (int i = r1, t = 0; i >= l1; i--, t++)
        {
            for (int j = r - t, k = 2; j >= l; j--, k *= 2)
            {
                if ((c >> j) & 1)
                {
                    ld x = 2 * pi; x /= k;
                    CPh(x, cbit, i);
                }
            }
        }
    }

    void CCFSUM(int c, int c1, int c2, int l1, int r1)
    {
        int l = 0;
        int r = r1 - l1;
        for (int i = r1, t = 0; i >= l1; i--, t++)
        {
            for (int j = r - t, k = 2; j >= l; j--, k *= 2)
            {
                if ((c >> j) & 1)
                {
                    ld x = 2 * pi; x /= k;
                    CCPh(x, c1, c2, i);
                }
            }
        }
    }

    void FSUB(int c, int l1, int r1)
    {
        int l = 0;
        int r = r1 - l1;
        for (int i = r1, t = 0; i >= l1; i--, t++)
        {
            for (int j = r - t, k = 2; j >= l; j--, k *= 2)
            {
                if ((c >> j) & 1)
                {
                    ld x = 2 * pi; x /= k;
                    Ph(-x, i);
                }
            }
        }
    }

    void CFSUB(int c, int cbit, int l1, int r1)
    {
        int l = 0;
        int r = r1 - l1;
        for (int i = r1, t = 0; i >= l1; i--, t++)
        {
            for (int j = r - t, k = 2; j >= l; j--, k *= 2)
            {
                if ((c >> j) & 1)
                {
                    ld x = 2 * pi; x /= k;
                    CPh(-x, cbit, i);
                }
            }
        }
    }

    void CCFSUB(int c, int c1, int c2, int l1, int r1)
    {
        int l = 0;
        int r = r1 - l1;
        for (int i = r1, t = 0; i >= l1; i--, t++)
        {
            for (int j = r - t, k = 2; j >= l; j--, k *= 2)
            {
                if ((c >> j) & 1)
                {
                    ld x = 2 * pi; x /= k;
                    CCPh(-x, c1, c2, i);
                }
            }
        }
    }

    void FSUM_mod_N(int a, int N, int c1, int c2, int l, int r, int zero_bit)
    {
        CCFSUM(a, c1, c2, l, r);
        FSUB(N, l, r);
        RQFT(l, r);
        CNOT(r, zero_bit);
        QFT(l, r);
        CFSUM(N, zero_bit, l, r);
        CCFSUB(a, c1, c2, l, r);
        RQFT(l, r);
        X(r);
        CNOT(r, zero_bit);
        X(r);
        QFT(l, r);
        CCFSUM(a, c1, c2, l, r);
    }

    void FSUB_mod_N(int a, int N, int c1, int c2, int l, int r, int zero_bit)
    {
        CCFSUB(a, c1, c2, l, r);
        RQFT(l, r);
        X(r);
        CNOT(r, zero_bit);
        X(r);
        QFT(l, r);
        CCFSUM(a, c1, c2, l, r);
        CFSUB(N, zero_bit, l, r);
        RQFT(l, r);
        CNOT(r, zero_bit);
        QFT(l, r);
        FSUM(N, l, r);
        CCFSUB(a, c1, c2, l, r);
    }


    void CMULT_mod_N(int a, int N, int c, int lx, int rx, int lb, int rb, int zero_bit)
    {
        QFT(lb, rb);
        for (int i = lx, mul = a % N; i <= rx; i++, mul = (((mul << 1) % N + N) % N))
        {
            FSUM_mod_N(mul, N, c, i, lb, rb, zero_bit);
        }
        RQFT(lb, rb);
    }

    void REV_CMULT_mod_N(int a, int N, int c, int lx, int rx, int lb, int rb, int zero_bit)
    {
        int n = rx - lx;
        QFT(lb, rb);
        for (int i = rx; i >= lx; i--, n--)
        {
            //(pow2(2, n, N) * a) % N
            FSUB_mod_N((pow2(2, n, N) * a) % N, N, c, i, lb, rb, zero_bit);
        }
        RQFT(lb, rb);
    }

    void SWAP(int c1, int c2)
    {
        CNOT(c1, c2);
        CNOT(c2, c1);
        CNOT(c1, c2);
    }

    void RANGE_SWAP(int l, int l1, int size)
    {
        for (int i = 0; i < size; ++i)
        {
            SWAP(l + i, l1 + i);
        }
    }

    void CSWAP(int control, int c1, int c2)
    {
        CCNOT(control, c1, c2);
        CCNOT(control, c2, c1);
        CCNOT(control, c1, c2);
    }

    void RANGE_CSWAP(int control, int l, int l1, int size)
    {
        for (int i = 0; i < size; ++i)
        {
            CSWAP(control, l + i, l1 + i);
        }
    }


    void CU(int a, int N, int c, int lx, int rx, int lb, int rb, int zero_bit)
    {
        CMULT_mod_N(a, N, c, lx, rx, lb, rb, zero_bit);
        RANGE_CSWAP(c, lx, lb, rx - lx + 1);
        REV_CMULT_mod_N(RevModNum(a, N), N, c, lx, rx, lb, rb, zero_bit);
    }

    //rc не включительно
    // QSHOR считает, что единица в состоянии уже поставлена
    int QSHOR(int a, int N, int rc, int lx, int rx, int lb, int rb, int zero_bit)
    {
        int ans = 0;
        for (int i = 0; i < rc; ++i)
            H(i);
        for (int i = 0; i < rc; i++)
        {
            CU(a, N, i, lx, rx, lb, rb, zero_bit);
            cout << "*";
            a = pow2(a, pow2(2, i, 1e9), N);
        }
        cout << "\n";
        RQFT(0, rc - 1);
        for (int i = 0; i < rc; ++i)
            Mes(i);
        for (int i = 0; i < 1 << size_in_bits; ++i)
            if (v[i].SqMod() > eps)
                ans = gcd(ans, (i >> rc) - 1);
        return ans;
    }

    void fast_QSHOR(int a, int N, int n, int lx, int rx, int lb, int rb, int zero_bit)
    {
        int buff = 0;
        int ans = 0;
        n <<= 1;

        // --- ZERO STEP ---
        H(0);
        CU(a, N, 0, lx, rx, lb, rb, zero_bit);
        H(0);
        int M = Mes(0);
        buff = M;
        if (M)
            X(0);
        // ---    END   ---

        for (int i = 1; i < n - 1; i++)
        {
            H(0);
            int value = pow2(a, pow2(2, i), N);
            CU(value, N, 0, lx, rx, lb, rb, zero_bit);
            for (int j = i - 1, k = 2; j >= 0; j--, k++)
            {
                if ((buff >> j) & 1)
                    R(k, 0);
            }
            M = Mes(0);
            buff |= M << i;
            if (M)
                X(0);
        }


        // --- LAST STEP ---
        H(0);
        int value = pow2(a, pow2(2, n - 1), N);
        CU(value, N, 0, lx, rx, lb, rb, zero_bit);
        for (int j = n - 2, k = 2; j >= 0; j--, k++)
        {
            if ((buff >> j) & 1)
                R(k, 0);
        }
        Mes(0);
        // ---    END   ---
    }


    void OutReal()
    {
        for (int i = 0; i < (1 << size_in_bits); ++i)
        {
            if (v[i].SqMod() > eps)
                cout << toBits(i, size_in_bits) << ": " << v[i] << std::endl;
        }
    }

    //изменяет состояние кубита
    int Mes(int ind)
    {
        GateCount++;
        ld p0, p1; p1 = p0 = 0;
        int flag = 0;
        std::random_device r;
        ld lower_bound = 0;
        ld upper_bound = 1;
        std::uniform_real_distribution<ld> unif(lower_bound, upper_bound);
        std::default_random_engine re(r());
        ld alpha = unif(re);

        for (int i = 0; i < v.size(); ++i)
        {
            if ((i & (1 << ind)) == 0)
            {
                p0 += v[i].SqMod();
            }
            else
            {
                p1 += v[i].SqMod();
            }
        }
        if (alpha <= p0)
            flag = 0;
        else
            flag = 1;

        ld x = flag ? p1 : p0;
        for (int i = 0; i < v.size(); ++i)
        {
            if (((i >> ind) & 1) == flag)
            {
                v[i] = v[i] / sqrt(x);
            }
            else
            {
                v[i] = complex(0, 0);
            }
        }

        return flag;
    }

    //выводит целочисленное представление двоичного состояния
    friend int MesAll(qMatrix& Q)
    {
        std::random_device r;
        ld lower_bound = 0;
        ld upper_bound = 1;
        std::uniform_real_distribution<ld> unif(lower_bound, upper_bound);
        std::default_random_engine re(r());
        ld alpha = unif(re);
        int ans = 0;

        for (int i = 0; i < Q.v.size(); ++i)
        {
            if (alpha <= Q[i].SqMod())
            {
                ans = i;
                break;
            }
            else
            {
                alpha -= Q[i].SqMod();
            }
        }

        for (int i = 0; i < Q.v.size(); ++i)
        {
            complex one(1, 0);
            complex zero(0, 0);
            if (i == ans)
                Q[i] = one;
            else
                Q[i] = zero;
        }

        return ans;
    }

    friend ostream& operator<<(ostream& out, qMatrix A)
    {
        int size = (1 << A.size_in_bits);
        for (int i = 0; i < size; ++i)
        {
            cout << toBits(i, A.size_in_bits) << ": " << A[i] << "\n";
        }
        return out;
    }

    friend void CARRY(qMatrix& q, int a1, int a2, int a3, int a4);
    friend void RCARRY(qMatrix& q, int a1, int a2, int a3, int a4);
    friend void SUM(qMatrix& q, int a1, int a2, int a3);

    friend void CARRY2(qMatrix& q, int a1, int a2, int a3, int a4, int a);
    friend void RCARRY2(qMatrix& q, int a1, int a2, int a3, int a4, int a);
    friend void SUM2(qMatrix& q, int a1, int a2, int a3, int a);

    friend void f1(qMatrix& q);
    friend void f2(qMatrix& q);
    friend void f3(qMatrix& q);
    friend void f4(qMatrix& q);
};
