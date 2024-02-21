#include "FuncTaylor.h"

#include <cmath>
#include <cstring>
using namespace std;

// constructor that transfer information from fa_ to data member terms.
FuncTaylor::FuncTaylor(const double *fa_, double a_, int order) : a(a_), terms(++order)
{
   if (order < 1)
   {
      cout << "Error in FuncTaylor constructor : order is inappropriate" << endl;
      cout << "Order parameter should be greater than -1" << endl;
      throw order;
   }
   memcpy((void *)&(terms(0)), (void *)fa_, sizeof(double) * order);
}

// basic printing function
void FuncTaylor::print(void) const
{
   const double zero = 0.;
   const int size = terms.getSize();
   if (size)
   {
      cout << terms(0);
      for (int i = 1; i < size; ++i)
      {
         cout << " , " << terms(i);
      }
   }
   cout << endl;
   return;
}
// get a data member
double FuncTaylor::getPoint(void) const { return a; }

// evaluate function for a given x
double FuncTaylor::operator()(const double &x) const
{
   double func_x = 0.;
   double t = x - a;
   for (int i = terms.getSize() - 1; i > 0; --i)
   {
      func_x += terms(i);  // add f^i(a) to the sum
      func_x *= t / i;     // multiply by (x-a)/i the sum
                           //      => rise by 1 order previous (x-a) term of the sum
                           //      => complete one more multiplication of factorial computation in divisor term of the sum
   }
   func_x += terms(0);  // last terms need to be treated apart to avoid division by 0

   return func_x;
}
// add two functions
FuncTaylor FuncTaylor::operator+(const FuncTaylor &p) const
{
   // only arithmetic of Taylor series around same points is possible in this class
   if (a != p.a)
   {
      cout << "+ operator works only in between Taylor expansion around the same point" << endl;
      throw -1;
   }
   const int size = terms.getSize();
   const int psize = p.terms.getSize();
   // checking size is mandatory to use properly Vector class.
   if (size > psize)
   {
      Vector terms_tmp = terms + p.terms;
      return FuncTaylor(&(terms_tmp(0)), a, size - 1);
   }
   else
   {
      Vector terms_tmp = p.terms + terms;
      return FuncTaylor(&(terms_tmp(0)), a, psize - 1);
   }
}
FuncTaylor FuncTaylor::operator-(const FuncTaylor &p) const
{
   // only arithmetic of Taylor series around same points is possible in this class
   if (a != p.a)
   {
      cout << "- operator works only in between Taylor expansion around the same point" << endl;
      throw -6;
   }
   const int size = terms.getSize();
   const int psize = p.terms.getSize();
   // checking size is mandatory to use properly Vector class.
   if (size < psize)
   {
      Vector terms_neg = p.terms - terms;
      // Use of *= is a little dirty here. Maybe in this case using directly a loop on terms may be more efficient.
      // In term of instruction number it is cheaper.
      terms_neg *= -1.;
      return FuncTaylor(&(terms_neg(0)), a, psize - 1);
   }
   else
   {
      Vector terms_neg = terms - p.terms;
      return FuncTaylor(&(terms_neg(0)), a, size - 1);
   }
}
FuncTaylor FuncTaylor::operator*(const FuncTaylor &p) const
{
   // only arithmetic of Taylor series around same points is possible in this class
   if (a != p.a)
   {
      cout << "* operator works only in between Taylor expansion around the same point" << endl;
      throw -6;
   }
   const double zero = 0.;
   const int size = terms.getSize();
   const int psize = p.terms.getSize();
   // Here we use "zero" local variable as an array of size 1. It avoid
   // creation of any extra specific array to create null Taylor expansion.
   if (size == 1 && terms(0) == zero) return FuncTaylor(&zero, a, 0);
   if (psize == 1 && p.terms(0) == zero) return FuncTaylor(&zero, a, 0);
   // compute product order and create a Vector instance for it with null coefficients.
   const int prod_order = size + psize - 2;
   const int prod_size = prod_order + 1;
   Vector terms_prod(prod_size);
   for (int i = 0; i < prod_size; ++i)
      terms_prod(i) = zero;  // terms_prod*=zero; may be a valuable solution (we exchange call to () operator by
                             //  multiplication by zero).

   // loop on p terms
   // only if p_terms_i is not null, we do the product with (*this). It cost 'psize' test
   // but it save 'size' multiplications and additions by zero when null.
   // initialize with constant term of p as there is no factorial correction in this case
   double p_terms_i = p.terms(0);
   if (p_terms_i != zero)
      for (int j = 0; j < size; ++j) terms_prod(j) += terms(j) * p_terms_i;
   for (int i = 1; i < psize; ++i)
   {
      p_terms_i = p.terms(i);
      if (p_terms_i != zero)
      {
         terms_prod(i) += terms(0) * p_terms_i;
         double fi = 1.;
         for (int j = 1; j < size; ++j)
         {
            // factorial impact:
            //    as a term at (i+j) index is implicitly supposed to be multiplied by (x-a)^(i+j)/(i+j)!
            //    the product at this index must be corrected so that this implicit convention is respected
            //    in the result
            //    terms(j)*(x-a)^j/j!  * p.terms(i)*(x-a)^i/i!
            //    = terms(j)*p.terms(i)*(x-a)^(i+j)/(i!*j!)
            //    = terms(j)*p.terms(i)*(x-a)^(i+j)/(i+j)! * (i+j)!/(i!*j!)
            //    = terms(j)*p.terms(i)*(x-a)^(i+j)/(i+j)! * (i+j)*(i+j-1)*...(i+1)/j!
            //
            //    Thus terms_prod(i+j) is augmented by terms(j)*p.terms(i)*(i+j)*(i+j-1)*...(i+1)/(j*(j-1)*...*1)
            //    which can be split as follows 
            //    terms(j)*p.terms(i)* (i+j)/j *(i+j-1)/(j-1) *...* (i+3)/3 * (i+2)/2 * (i+1)/1 
            //    or
            //    terms(j)*p.terms(i)* (i/j+1) * (i/(j-1)+1) *...* (i/3+1) * (i/2+1) * (i+1) 
            //
            fi *= (1.+(1.*i)/(1.*j));
            terms_prod(i + j) += terms(j) * p_terms_i * fi;
         }
      }
   }

   return FuncTaylor(&(terms_prod(0)), a, prod_order);
}

FuncTaylor FuncTaylor::deriv(const unsigned int &nth) const
{
   int size = terms.getSize();
   int nsize = size - nth;
   nsize = (nsize > 0) ? nsize : 1;
   Vector terms_derived(nsize);
   terms_derived(0)=0.;
   for (int i = nth; i < size; ++i) terms_derived(i - nth) = terms(i);
   return FuncTaylor(&(terms_derived(0)), a, nsize-1);
}

// safety
FuncTaylor &FuncTaylor::operator=(const FuncTaylor &p)
{
   cout << "As this class is just for an exercise without changing Vector class this = operator is not implementable correctly"
        << endl;
   cout
       << "All it can do by default is using = operator of Vector class which is quite restrictive due to the left hand size rule"
       << endl;
   cout << "Don't use it" << endl;
   throw -2;
}
