#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
// mandatory inclusion of class declaration
#include "Real.h"

using namespace std;

// A example of macro usage (rather complicate but
// give macro EXPP and EXPM that give double constant
// syntax as function of exponent)
#define MX 15
#define MXM1 14
#define AFTERP(a) 1.e ## a
#define AFTERM(a) 1.e- ## a
#define EXPP(a)  AFTERP(a)
#define EXPM(a)  AFTERM(a)

// Constructor
Real::Real(double m_, int e_) : m(m_),e(e_)
{
    const double one = 1.;
    const double mone = -1.;
    const double zone = 0.1;
    const double mzone = -0.1;
    const double eps=EXPM(MX);
    if (m < eps && m>-eps)
    {
        m=0.;
        e = 0;
    }
    else
    {
        while ( !( m < one ) || !( m > mone ) )
        {
            m /= 10.;
            ++e;
        }
        while ( m < zone  && m > mzone)
        {
            m=floor(m*EXPP(MX))/EXPP(MXM1); //equivalent to "m *= 10.;" but truncate to MX digit so that approximation on digit 16 and higher do 
                                            // not interfere in computation. This is improving a little precision of the class. 
            --e;
        }
    }
}

// display function
void Real::display() const
{
    cout <<setprecision(MX)<< m << "E";
    cout << e << endl;
}

// addition operator
Real Real::operator+(const Real & rhs) const
{
    const double zero=0.;
    // we must check if one operand is null because a zero value should not
    // hide a tiny one. Only need to test mantissa with zero as constructor already
    // do the job to compare with EXPM(MX)
    // Test first "this"
    if (m == zero)
    {
        return rhs;
    }
    // Test second rhs
    if (rhs.m == zero)
    {
        return *this;
    }
    // Now  use exponent
    int de = e-rhs.e;
    if (abs(de) > MX)
    {
        if (de > 0) return *this;
        else return rhs;
    }
    else
    {
        if (de > 0)
        {
            return Real(m*pow(10.,de)+rhs.m,rhs.e);
        }
        else
        {
            return Real(m+rhs.m*pow(10.,-de),e);
        }

    }
}

// subtraction operator
Real Real::operator-(const Real & rhs) const
{
    const double zero=0.;
    // we must check if one operand is null because a zero value should not
    // hide a tiny one. Only need to test mantissa with zero as constructor already
    // do the job to compare with EXPM(MX)
    // Test first "this"
    if (m == zero)
    {
        return Real(-rhs.m, rhs.e);
    }
    // Test second rhs
    if (rhs.m == zero)
    {
        return *this;
    }
    // Now  use exponent
    int de = e-rhs.e;
    if (abs(de) > MX)
    {
        if (de > 0) return *this;
        else return Real(-rhs.m, rhs.e);
    }
    else
    {
        if (de > 0)
        {
            return Real(m*pow(10.,de)-rhs.m,rhs.e);
        }
        else
        {
            return Real(m-rhs.m*pow(10.,-de),e);
        }

    }
}

// multiplication operator
Real Real::operator*(const Real & rhs) const
{
   // Test to check that e+rh.e is not overpassing the numeric representation of an int.
   if (e > 0 && rhs.e > 0 && (e + rhs.e) < 0)
   {
      cout << "error : e+rhs.e is to big: " << e << "+" << rhs.e
           <<" greater then int range " << endl;
          throw -23;
   }
   return Real(m * rhs.m, e + rhs.e);
}

// division operator
Real Real::operator/(const Real & rhs) const
{
    if (rhs.m == 0.)
    {
        cout << "Warning : trying to divide by a null real number !!! " << endl;
        throw -1;
    }
    else
    {
        return Real( m/rhs.m,e-rhs.e);
    }
}

// power of n
Real Real::power(int n) const
{
   Real prod(*this);
   Real fact(*this);

   // if negative power we take 1/(*this) for the product
   if (n < 0)
   {
      fact = Real(1., 0) / (*this);
      prod = fact;
      n = -n;
    }
    // if null power we skip the product and give directely the result : 1
    else if (n < 1)
    {
       prod = Real(1., 0);
       n=1;
    }

    // product to replace power
    // basic implementation but a log n implementation exist
    while (--n)
    {
       prod = prod * fact;
    }
    return prod;
}

// exponetial
Real Real::exp() const
{
   // exp^(m.10^e)=(exp^(10^e))^m=(k.10^l)^m=k^m.(10^m)^l
   // with k.10^l found by computing (..(((exp^10)^10)^10)....) e times if e>0
   //               or  by computing (..(((exp^(1/10))^(1/10))^(1/10))....) |e| times if e<0
   //               or by taking k=1,l=0 if e=0 => result is just exp^m
   Real ex(std::exp(10.), 0);
   if (e < 0)
   {
      ex = Real(std::exp(0.1), 0);
      for (int i = 1,n=-e; i < n; ++i)
         ex = Real(pow(ex.m, 0.1) * pow(pow(10., ex.e), 0.1), 0);
   }
   else if (e < 1)
      return Real(std::exp(m), 0);
   else
      for (int i = 1; i < e; ++i) ex = ex.power(10);
   Real dm(pow(10., m));
   Real km(pow(ex.m, m));
   return km * (dm.power(ex.e));
}
