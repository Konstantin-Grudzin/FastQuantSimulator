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
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>
#include<omp.h>

using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;

using namespace std;

typedef long double ld;
static ld aa = 1 / sqrt(2);
static ld pi = acos(-1);
static ld eps = 1e-5;

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

class Complex
{
	ld a;
	ld b;
public:
	Complex(ld _a = 0, ld _b = 0) : a(_a), b(_b)
	{

	}
	Complex(complex<double> y)
	{
		a = y.real();
		b = y.imag();
	}
	Complex operator+(Complex y)
	{
		Complex z(a + y.a, b + y.b);
		return z;
	}
	Complex operator-(Complex y)
	{
		Complex z(a - y.a, b - y.b);
		return z;
	}
	Complex operator*(Complex y)
	{
		Complex z(a * y.a - b * y.b, a * y.b + b * y.a);
		return z;
	}
	Complex operator/ (ld x)
	{
		Complex z(a / x, b / x);
		return z;
	}
	Complex operator* (ld x)
	{
		Complex z(a * x, b * x);
		return z;
	}
	Complex operator* (complex<double> y)
	{
		Complex z(0, 0);
		z = y;
		return (*this) * z;
	}
	Complex& operator*= (Complex y)
	{
		Complex z(a * y.a - b * y.b, a * y.b + b * y.a);
		this -> operator=(z);
		return *this;
	}
	Complex operator= (std::complex<double> cmp)
	{
		a = cmp.real();
		b = cmp.imag();
		return *this;
	}
	Complex& operator+=(Complex y)
	{
		a += y.a;
		b += y.b;
		return *this;
	}
	bool operator==(Complex v)
	{
		return (a == v.a) && (b == v.b);
	}
	bool operator!=(Complex v)
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

	friend istream& operator>>(istream& in, Complex a)
	{
		in >> a.a >> a.b;
		return in;
	}

	friend ostream& operator<<(ostream& out, Complex a)
	{
		out << a.a;
		if (abs(a.b) <= eps)
		{
			return out;
		}

		if (a.b > 0)
			out << " + " << a.b << " * i";
		else
			out << " - " << abs(a.b) << " * i";
		return out;
	}

};

static Complex exp_i_phi(double phi)
{
	Complex a(cos(phi), sin(phi));
	return a;
}

static Complex obrk2(aa);

//всегда начинаем из 0000
class qMatrix
{
	int size_in_bits;
	atomic<int> GateCount = 0;
	vector<Complex> v;
	vector<Complex> buffer;
	Complex* bufptr, * vptr;
	ld error;
	bool random_bit;
public:
	qMatrix(int size = 1, ld err = 0.00000000000001, bool rb = false) : size_in_bits(size), error(err), random_bit(rb)
	{
		v.resize(1 << (size_in_bits));
		buffer.resize(1 << (size_in_bits));
		bufptr = buffer.data();
		vptr = v.data();
		v[0] = Complex(1);
	}

	inline int GetGateCount()
	{
		return GateCount;
	}

	Complex& operator[](int ind)
	{
		return v[ind];
	}

	void reset()
	{
		GateCount = 0;
		v[0] = Complex(1);
		auto* tmp = v.data();
#pragma omp parallel for
		for (int i = 1; i < v.size(); ++i)
		{
			tmp[i] = Complex(0);
		}
	}

	const qMatrix& operator=(const vector<Complex>& v1)
	{
		if (v.size() != v1.size())
		{
			throw std::runtime_error("Enequal size!");
		}
		auto* ptr = v.data();
		auto* ptr1 = v1.data();
#pragma omp parallel for
		for (int i = 0; i < v1.size(); ++i)
			ptr[i] = ptr1[i];
		return *this;
	}

	// initialize vector state b in qMatrix int [start,finish) qubits
	void initialize(Vector& b, int start, int finish)
	{
		if (pow(2, finish - start) != b.rows())
			throw "Wrong init!";
		for (int i = start; i < b.rows(); ++i)
		{
			v[i] = b(i) * b(i);
		}
	}

	// initialize vector state b in qMatrix int [0,finish) qubits
	void initialize(Vector& b, int finish)
	{
		initialize(b, 0, finish);
	}

	// Use UnitaryGate on [start,finish) qubits
	void UnitaryGate(Matrix& U, int start, int finish)
	{
		for (int i = 0; i < (1 << size_in_bits); ++i)
			bufptr[i] = vptr[i];
		int mask = (1 << finish) - 1;
		mask >>= start;
		mask <<= start;
		int rmask = ~mask;
		//#pragma omp parallel for
		for (int i = 0; i < v.size(); ++i)
		{
			int pos = i & mask;
			int nind = i & rmask;
			for (int j = 0; j < U.cols(); ++j)
			{
				vptr[i] += bufptr[nind + (j << start)] * U(pos, j);
			}
		}
	}

	// Use UnitaryGate on [0,finish) qubits
	void UnitaryGate(Matrix& U, int finish)
	{
		UnitaryGate(U, 0, finish);
	}

	// Add Control to Ugate [start,finish)
	void CUnitaryGate(Matrix& U, int cbit, int start, int finish)
	{
		for (int i = 0; i < (1 << size_in_bits); ++i)
			bufptr[i] = vptr[i];
		int mask = (1 << finish) - 1;
		mask >>= start;
		mask <<= start;
		int rmask = ~mask;
		//#pragma omp parallel for
		for (int i = 0; i < v.size(); ++i)
		{
			if (i & (1 << cbit))
			{
				int pos = i & mask;
				int nind = i & rmask;
				for (int j = 0; j < U.cols(); ++j)
				{
					vptr[i] += bufptr[nind + (j << start)] * U(pos, j);
				}
			}
		}
	}

	// Add Control to Ugate [0,finish)
	void CUnitaryGate(Matrix& U, int cbit, int finish)
	{
		CUnitaryGate(U, cbit, 0, finish);
	}

	//Make QPE use control qubits [start1,finish1) to UnitaryGate on qubits [start,finsh)
	void QPE(Matrix U, int start, int finish, int start1, int finish1)
	{
		// construct QPE
		for (int i = start1; i < finish1; ++i)
		{
			H(i);
		}
		for (int i = start1; i < finish1; ++i,U *= U)
		{
			CUnitaryGate(U, i, start,finish);
		}
		RQFT(start1, finish1 - 1);
	}

	void U3(ld theta, ld phi, ld lamda, int ind)
	{
		Complex cf1 = cos(theta / 2), cf2 = exp_i_phi(lamda) * sin(-theta / 2);
		Complex cf3 = exp_i_phi(phi) * sin(theta / 2), cf4 = exp_i_phi(phi + lamda) * cos(theta / 2);
#pragma omp parallel for
		for (int i = 0; i < v.size(); ++i)
		{
			bufptr[i] = vptr[i];
		}
#pragma omp parallel for
		for (int i = 0; i < v.size(); ++i)
		{
			if (i & (1 << ind))
				vptr[i] = cf3 * bufptr[i - (1 << ind)] + cf4 * bufptr[i];
			else
				vptr[i] = cf1 * bufptr[i] + cf2 * bufptr[i + (1 << ind)];
		}
	}

	void CU3(ld theta, ld phi, ld lamda, int cbit, int ind)
	{
		Complex cf1 = cos(theta / 2), cf2 = exp_i_phi(lamda) * sin(-theta / 2);
		Complex cf3 = exp_i_phi(phi) * sin(theta / 2), cf4 = exp_i_phi(phi + lamda) * cos(theta / 2);
#pragma omp parallel for
		for (int i = 0; i < v.size(); ++i)
		{
			bufptr[i] = vptr[i];
		}
#pragma omp parallel for
		for (int i = 0; i < v.size(); ++i)
		{
			if (i & (1 << cbit))
			{
				if (i & (1 << ind))
					vptr[i] = cf3 * bufptr[i - (1 << ind)] + cf4 * bufptr[i];
				else
					vptr[i] = cf1 * bufptr[i] + cf2 * bufptr[i + (1 << ind)];
			}
		}
	}

	void Rx(ld theta, int ind)
	{
		U3(theta, -pi / 2, pi / 2, ind);
	}

	void Ry(ld theta, int ind)
	{
		U3(theta, 0, 0, ind);
	}

	void Rz(ld lamda, int ind)
	{
		Complex cf1 = exp_i_phi(-lamda / 2), cf2 = exp_i_phi(lamda / 2);
#pragma omp parallel for
		for (int i = 0; i < v.size(); ++i)
		{
			if (i & (1 << ind))
				vptr[i] *= cf2;
			else
				vptr[i] *= cf1;
		}
	}

	void X(int ind)
	{
		GateCount++;
		int size = (1 << size_in_bits);
#pragma omp parallel for
		for (int i = 0; i < size; ++i)
		{
			if (i & (1 << ind))
			{
				swap(v[i], v[i ^ 1 << ind]);
			}
		}
	}

	void Z(int ind)
	{
		GateCount++;
		auto* ptr = v.data();
#pragma omp parallel for
		for (int i = 0; i < (1 << size_in_bits); ++i)
		{
			if (i & (1 << ind))
				ptr[i] = ptr[i] * (-1);
		}
	}

	void H(int ind)
	{
#pragma omp parallel for
		for (int i = 0; i < 1 << (1 << (size_in_bits - ind - 1)); ++i)
		{
			for (int j = i << (ind + 1); j < (1 << ind)+(i<<(ind+1)); ++j)
			{
				auto t = vptr[i] * obrk2;
				auto t1 = vptr[i + (1 << ind)] * obrk2;
				vptr[i] = t - t1;
				vptr[i + (1 << ind)] = t + t1;
			}
		}
		//for (int i = 0; i < v.size(); ++i)
		//{
		//	if (i & (1 << ind))
		//	{
		//		auto t = vptr[i]*obrk2;
		//		auto t1 = vptr[i - (1 << ind)]*obrk2;
		//		vptr[i] = t1-t;
		//		vptr[i + (1 << ind)] = t + t1;
		//	}
		//}
	}

	void CNOT(int cind, int chnind)
	{
		//TODO: extractin' for by 1<<(n-2)
		GateCount++;
		qMatrix a(size_in_bits);
		int size = 1 << (size_in_bits);
		for (int i = 0; i < size; ++i)
		{
			if (i & (1 << cind))
			{
				buffer[i ^ (1 << chnind)] = v[i];
			}
			else
			{
				buffer[i] = v[i];
			}
		}
		*this = buffer;
	}

	void CCNOT(int cind0, int cind1, int chnind)
	{
		GateCount++;
		int size = 1 << (size_in_bits);
		for (int i = 0; i < size; ++i)
		{
			if ((i & (1 << cind0)) && (i & (1 << cind1)))
			{
				buffer[i ^ (1 << chnind)] = v[i];
			}
			else
			{
				buffer[i] = v[i];
			}
		}
		*this = buffer;
	}

	void Ph(double theta, int bit)
	{
		GateCount++;
		ld alpha;
		int size = 1 << (size_in_bits);
		if (random_bit)
		{
			std::random_device r;
			ld lower_bound = theta * (1 - error);
			ld upper_bound = theta * (1 + error);
			std::uniform_real_distribution<ld> unif(lower_bound, upper_bound);
			std::default_random_engine re(r());
			alpha = unif(re);
		}
		else
		{
			alpha = theta;
		}
		Complex exp = exp_i_phi(alpha);
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
		Complex exp = exp_i_phi(theta);
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
		Complex exp = exp_i_phi(theta);
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

	int fast_QSHOR(int a, int N, int n, int lx, int rx, int lb, int rb, int zero_bit)
	{
		int buff = 0;
		n <<= 1;

		// --- PREPARATION ---
		X(1);
		// -------------------

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
			//cout << "*";
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
		int ans = 0;
		for (int i = 0; i < 1 << size_in_bits; ++i)
			if (v[i].SqMod() > eps)
			{
				ans = gcd(ans, (i >> 1) - 1);
				cout << (i >> 1) - 1 << "\n";
			}
		return ans;
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
#pragma omp parallel for
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
#pragma omp parallel for
		for (int i = 0; i < v.size(); ++i)
		{
			if (((i >> ind) & 1) == flag)
			{
				v[i] = v[i] / sqrt(x);
			}
			else
			{
				v[i] = Complex(0);
			}
		}

		return flag;
	}

	//If you need to measure some special value on qubit 
	//for some reason
	int CheatMes(int ind, int value)
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
#pragma omp parallel for
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
		if (value == 0)
			flag = 0;
		else
			flag = 1;
		ld x = flag ? p1 : p0;
#pragma omp parallel for
		for (int i = 0; i < v.size(); ++i)
		{
			if (((i >> ind) & 1) == flag)
			{
				v[i] = v[i] / sqrt(x);
			}
			else
			{
				v[i] = Complex(0);
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
			Complex one(1, 0);
			Complex zero(0, 0);
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