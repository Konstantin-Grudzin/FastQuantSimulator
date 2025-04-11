#include "qMatrix.h"

void f1(qMatrix& q)
{

}
void f2(qMatrix& q)
{
    q.X(1);
}
void f3(qMatrix& q)
{
    q.CNOT(0, 1);
}
void f4(qMatrix& q)
{
    q.X(0);
    q.CNOT(0, 1);
    q.X(0);
}


void CARRY(qMatrix& q, int a1, int a2, int a3, int a4)
{
    q.CCNOT(a2, a3, a4);
    q.CNOT(a2, a3);
    q.CCNOT(a1, a3, a4);
}
void RCARRY(qMatrix& q, int a1, int a2, int a3, int a4)
{
    q.CCNOT(a1, a3, a4);
    q.CNOT(a2, a3);
    q.CCNOT(a2, a3, a4);
}
void SUM(qMatrix& q, int a1, int a2, int a3)
{
    q.CNOT(a2, a3);
    q.CNOT(a1, a3);
}



void CARRY2(qMatrix& q, int a1, int a2, int a3, int a4, int a)
{
    a1 -= 3;
    a3 -= 3;
    a4 -= 3;
    if ((a >> a2) & 1)
        q.CNOT(a3, a4);
    if ((a >> a2) & 1)
        q.X(a3);
    q.CCNOT(a1, a3, a4);
}
void RCARRY2(qMatrix& q, int a1, int a2, int a3, int a4, int a)
{
    a1 -= 3;
    a3 -= 3;
    a4 -= 3;
    q.CCNOT(a1, a3, a4);
    if ((a >> a2) & 1)
        q.X(a3);
    if ((a >> a2) & 1)
        q.CNOT(a3, a4);
}
void SUM2(qMatrix& q, int a1, int a2, int a3, int a)
{
    a1 -= 3;
    a3 -= 3;
    if ((a >> a2) & 1)
        q.X(a3);
    q.CNOT(a1, a3);
}