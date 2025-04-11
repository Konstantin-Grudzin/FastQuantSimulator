#define _CRT_SECURE_NO_WARNINGS

#include "qMatrix.h"
#include<ctime>
#include <chrono>


using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;

int Get_Biggest_pow_of_2(int x)
{
	int cnt = 0;
	int ans = 1;
	while (ans < x)
		ans <<= 1, cnt++;
	return ans;
}

//Реализация с 3n кубитами
void Three_n_Qubit_Addition()
{

	qMatrix q(10);
	cout << "Input two nubers with space:";

	int a, b; cin >> a >> b;
	for (int i = 0; i < 3; ++i)
	{
		if ((a >> i) & 1)
		{
			q.X(i);
		}
	}
	for (int i = 0; i < 3; ++i)
	{
		if ((b >> i) & 1)
		{
			q.X(6 + i);
		}
	}
	cout << "Thats how it's looks like in quants: "; q.OutReal();
	CARRY(q, 3, 0, 6, 4);
	CARRY(q, 4, 1, 7, 5);
	CARRY(q, 5, 2, 8, 9);
	q.CNOT(2, 8);
	SUM(q, 5, 2, 8);
	RCARRY(q, 4, 1, 7, 5);
	SUM(q, 4, 1, 7);
	RCARRY(q, 3, 0, 6, 4);
	SUM(q, 3, 0, 6);

	cout << "Summary of our addition: ";
	q.OutReal();
}

//Реализация с 2n кубитами
void Two_n_Qubit_Addition()
{
	qMatrix q(7);
	cout << "Input two nubers with space:";

	int a, b; cin >> a >> b;
	for (int i = 0; i < 3; ++i)
	{
		if ((b >> i) & 1)
		{
			q.X(3 + i);
		}
	}

	cout << "Thats how it's looks like in quants: "; q.OutReal();
	CARRY2(q, 3, 0, 6, 4, a);
	CARRY2(q, 4, 1, 7, 5, a);
	CARRY2(q, 5, 2, 8, 9, a);
	if ((a >> 2) & 1)
		q.X(5);
	SUM2(q, 5, 2, 8, a);
	RCARRY2(q, 4, 1, 7, 5, a);
	SUM2(q, 4, 1, 7, a);
	RCARRY2(q, 3, 0, 6, 4, a);
	SUM2(q, 3, 0, 6, a);
	cout << "Summary of our addition: ";
	q.OutReal();
}

void Deutsch(int func = 1)
{
	qMatrix q(2);
	vector<int> stat(4);
	cout << func << " is your function\n";
	for (int i = 0; i < 10000; ++i)
	{
		q.reset();
		q.H(0);
		//q.OutReal();cout << endl;
		q.X(1);
		//q.OutReal(); cout << endl;
		q.H(1);
		//q.OutReal(); cout << endl;

		if (func == 1)
			f1(q);
		else if (func == 2)
			f2(q);
		else if (func == 3)
			f3(q);
		else
			f4(q);

		q.H(0);
		//q.OutReal(); cout << endl;
		q.Mes(0);
		//q.OutReal(); cout << endl;
		stat[MesAll(q)]++;
	}
	for (int i = 0; i < 4; ++i)
	{
		cout << toBits(i, 2) << " " << stat[i] << endl;
	}
	cout << endl;
	bool bal = stat[1] && stat[3];
	bool sta = stat[0] && stat[2];
	if (bal)
	{
		cout << "Function is balanced\n";
	}
	else if (sta)
	{
		cout << "Function is static\n";
	}
	else
	{
		cout << "Neither\n";
	}
	cout << endl;
}

Matrix exp_iAt(Matrix A, double t = 150)
{
	std::complex<double> i(0, 1);
	Eigen::ComplexEigenSolver<Matrix> ces(A);
	Matrix V = ces.eigenvectors();
	Vector D = ces.eigenvalues();
	Matrix expA(D.size(), D.size());
	expA.setZero();
	for (int ind = 0; ind < D.size(); ++ind)
	{
		auto tmp = exp(i * D(ind) * t) * V.col(ind) * V.col(ind).transpose();
		expA += tmp;
	}
	return expA;
}

Vector GetEigenValues(const Matrix& A)
{
	Eigen::ComplexEigenSolver<Matrix> ces(A);
	Vector D = ces.eigenvalues();

	return D;
}

void HHL(int sz, Matrix& A, Vector& b)
{
	//set vars
	// C choose by wanted accuracy
	ld C = 1/1024, t = 150;
	int n_eig=8, n;
	//Create Ermit Matrix and Vector
	int power = Get_Biggest_pow_of_2(sz)*2;
	n = log2l(power);
	Matrix A_p(power, power);
	A_p.setZero();
	A_p.block(0, A.cols(), A.rows(), A.cols()) = A;
	A_p.block(A.rows(), 0, A.rows(), A.cols()) = A.conjugate().transpose();
	for (int i = 2 * sz; i < power; ++i)
	{
		A_p(i, i) = 1;
	}
	cout << A_p << endl;
	Vector l = GetEigenValues(A_p);
	double mmx = -numeric_limits<double>::infinity();
	for (int i = 0; i < l.rows(); ++i)
	{
		mmx = max(mmx, l(i).real());
	}
	A_p /= mmx;
	l /= mmx;
	//cout << A_p << endl;
	//eigenvalues
	Matrix A_exp = exp_iAt(A_p, t);
	Vector b_p(power);
	b_p.setZero();
	for (int i = 0; i < sz; ++i)
		b_p(i) = b(i);
	b_p.normalize();

	//First for measure (to make correspondence)
	qMatrix q(n + n_eig + 1);
	q.initialize(b_p, n);
	q.QPE(A_exp, 0, n, n, n + n_eig);
	q.OutReal();

	//auto top_ev_bin = get_top_ev_bins(q);
}

int main()
{
	//freopen("data.out", "w", stdout);
	ios_base::sync_with_stdio(0);
	cout.tie(NULL);
	random_device r;
	mt19937 x(r());
	Matrix A(2, 2);
	A << 1, 0.6, 0.6, 1;
	Vector b(2);
	b << 0, 1;
	const auto start{ std::chrono::steady_clock::now() };
	HHL(A.rows(), A, b);
	const auto end{ std::chrono::steady_clock::now() };
	const std::chrono::duration<double> elapsed_seconds{ end - start };
	cout << elapsed_seconds.count();
}
//working example
//cout << q.fast_QSHOR(11, 5 * 23, 8, 1, 8, 9, 17, 18) << "\n";