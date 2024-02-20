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
    //store a bool whether x is negative or not
    bool neg = (x < 0.);
    //take the absolute value of x
    double y = neg ? -x : x;
    //same computation as before taylor_expansion
    for (int i = 1; i < N; ++i) {
        a *= y / i;
        res += a;
    }
    //if x is negative, the result is the inverse of the result for the absolute value of x
    // => res=e^y=e^-x=(e^x)^-1=1./e^x
    // => exp(-x) = 1./exp(x)
    if (neg) res = 1./res;
    
    return res;
}

int main(int argc, char *argv[]) 
{
    //variable declaration
    double xmin, xmax, dx;
    int Order, nbpoints;
    // entering user input
    cout <<  "Give the range of x (xmin, xmax): ";
    cin >> xmin >> xmax;
    cout <<  "Give the number of points(not counting boundary points): "<<
    endl;
    cin >> nbpoints;
    nbpoints += 2;
    cout <<  "Give the order of the expansion: "<<endl;
    cin >> Order;

    //compute the step size
    dx = (xmax - xmin) / (nbpoints - 1);
    //loop over the range of x
    cout << "Taylor expansion of exp(x) from " << xmin << " to " << xmax << " with " << Order << " terms" << endl;
    for (int i = 0; i < nbpoints; ++i) {
        double x = xmin + i * dx;
        cout << "x " << i << " : " << taylor_expansion(x, Order) << endl;
    }
    cout << endl;

    //loop over the range of x to compute exact exponential value, the asked function, the recursive implementation and the negative implementation
    cout << "Comparison"<< endl;
    for (int i = 0; i < nbpoints; ++i) {
        double x = xmin + i * dx;
        cout << scientific << x << " : " << exp(x) << " : " << taylor_expansion(x, Order) << " : " << taylor_exp_r(x, Order) << " : " << taylor_exp_n(x, Order) << endl;
    }

    return 0;
}
