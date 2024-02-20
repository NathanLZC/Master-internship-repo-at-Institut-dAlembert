#include <iostream>
#include <cmath>
using namespace std;

double taylor_expansion(double x, int N) 
{
    double res = 1.;//for taylor expansion
    double a = 1.;//to compute x^k / k!
    //By increasing the value of N, the test below in the for statement is a simple inequality. This is cheaper
    //than the <= operator which is equivalent to the < operator and = test( = test itself is equivalent to not < and not >). 
    //Here N is not that big, becasue it represents the order of the taylor expansion.

    N++;

    for (int i = 1; i < N; ++i) {
        a *= x / i;
        res += a;
    }
    return res;
}

double recur_expapprox(double, int, double, double, int); 
//Here recursion is hidden from the user, and the user only needs to provide the first two arguments.
double taylor_exp_r(double x, int N) 
{
    if (N > 0) return recur_expapprox(x, 1, 1., 1., N+1); //this first call initiated the recursion
    else return 1.;
}

double recur_expapprox(double x, int i, double a, double res, int N) 
{
    //For the current order, compute the expansion
    a *= x/(i++);
    res += a;
    //Testing condition to stop the recursion
    //Here when i is equal to N, the recursion stops and the result is returned.
    if (i < N) return recur_expapprox(x, i, a, res, N);
    return res;
}

double taylor_exp_n(double x, int N) 
{
    double res = 1.;
    double a = 1.;
    N++;
    //
    bool neg = (x < 0.);
    //
    for (int i = 1; i < N; ++i) {
        a *= x / i;
        res += a;
    }
    return res;
}