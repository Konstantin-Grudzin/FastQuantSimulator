#define _CRT_SECURE_NO_WARNINGS

#include "qMatrix.h"
#include<ctime>


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

void Deutsch(int func=1)
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

int main()
{
    //freopen("data.out", "w", stdout);
    ios_base::sync_with_stdio(0);
    cout.tie(NULL);
    random_device r;
    mt19937 x(r());
    int func = x() % 4 + 1;
    Deutsch(1);
    Deutsch(2);
    Deutsch(3);
    Deutsch(4);
   /* Deutsch(func);*/
}