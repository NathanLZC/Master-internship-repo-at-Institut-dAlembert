#include <cmath>
#include <iostream>
using namespace std;

#include "FuncTaylor.h"

#define MXORDER 200

// compute derivative of polynomial function
// f(x)=-1+6x-2x^3+9x^4+3x^5
void deri_poly(double *f, const double &a, const int &n)
{
   for (int i = 0; i <= n; i++)
   {
      switch (i)
      {
         case 0:
         {
            f[i] = -1. + 6. * a - 2. * a * a * a + 9. * pow(a, 4) + 3. * pow(a, 5);
            break;
         }
         case 1:
         {
            f[i] = 6. - 6. * a * a + 36. * pow(a, 3) + 15. * pow(a, 4);
            break;
         }
         case 2:
         {
            f[i] = -12. * a + 108. * pow(a, 2) + 60. * pow(a, 3);
            break;
         }
         case 3:
         {
            f[i] = -12. + 216. * a + 180. * a * a;
            break;
         }
         case 4:
         {
            f[i] = 216. + 360. * a;
            break;
         }
         case 5:
         {
            f[i] = 360.;
            break;
         }
         default:
         {
            f[i] = 0.;
         }
      };
   }
}

int main()
{
   // A buffer to store f(a),f'(a),f''(a),....
   // Warning : in this program, we suppose that the maximum order for Taylor series is MXORDER
   double fv[MXORDER] = {0.};

   // forcing same output precision for everyone (hoping no one change this setting in their class)
   cout.precision(5);
   cout << scientific;

   // Instantiate 3 objects
   int n;         // order of highest term of Taylor expansion description
   double a;      // point at which Taylor expansion terms are computed
   double a_ref;  // point at which Taylor expansion terms are computed for all function to test arithmetic

   // Instantiate a object that will represents cos function :
   cout << "Input order of Taylor expansion to represent cosine function (less then " << MXORDER << ")" << endl;
   cin >> n;
   cout << "Input point at which Taylor expansion terms are computed (reference point)" << endl;
   cin >> a;
   a_ref = a;

   // compute derivative of cosine function
   for (int i = 0, k = 0; i <= n; i++)
   {
      double v;
      switch (k)
      {
         case 0:
         {
            v = cos(a);
            break;
         }
         case 1:
         {
            v = -sin(a);
            break;
         }
         case 2:
         {
            v = -cos(a);
            break;
         }
         case 3:
         {
            v = sin(a);
            break;
         }
      };
      k = (++k) % 4;
      fv[i] = v;
   }
   // create instance
   FuncTaylor f_cos(fv, a, n);

   // Instantiate a object that will represents exponential function :
   cout << "Input order of Taylor expansion to represent exponential function (less then " << MXORDER << ")" << endl;
   cin >> n;
   cout << "Input point at which Taylor expansion terms are computed" << endl;
   cin >> a;

   // compute derivative of exponential function
   double expa = exp(a);
   for (int i = 0; i <= n; i++) fv[i] = expa;  // all derivatives of exp function are exp function itself
   // create instance
   FuncTaylor f_exp(fv, a, n);
   if (a != a_ref)
   {
      expa = exp(a_ref);
      for (int i = 0; i <= n; i++) fv[i] = expa;
   }
   FuncTaylor f_exp_ref(fv, a_ref, n);

   // Instantiate a object that will represents 1/(1-x) function :
   cout << "Input order of Taylor expansion to represent 1/(1-x) function (less then " << MXORDER << ")" << endl;
   cin >> n;
   cout << "Input point at which Taylor expansion terms are computed" << endl;
   cin >> a;

   // compute derivative of 1/(1-x) function
   double c = 1. / (1. - a);
   fv[0] = c;
   for (int i = 1; i <= n; i++)
   {
      c *= i / (1. - a);
      fv[i] = c;  // derivatives are (1-a)^(-(n+1)).n!
   }

   // create instance
   FuncTaylor f_geom(fv, a, n);
   if (a != a_ref)
   {
      double c = 1. / (1. - a_ref);
      fv[0] = c;
      for (int i = 1; i <= n; i++)
      {
         c *= i / (1. - a_ref);
         fv[i] = c;
      }
   }
   FuncTaylor f_geom_ref(fv, a_ref, n);

   // Instantiate a object that will represents a polynomial function of order 5 :
   // f(x)=-1+6x-2x^3+9x^4+3x^5
   cout << "Input order of Taylor expansion to represent poly function (less then " << MXORDER << ")" << endl;
   cin >> n;
   cout << "Input point at which Taylor expansion terms are computed" << endl;
   cin >> a;

   // compute derivative of polynomial function
   deri_poly(fv, a, n);

   // create instance
   FuncTaylor f_poly(fv, a, n);
   if (a != a_ref) deri_poly(fv, a_ref, n);
   FuncTaylor f_poly_ref(fv, a_ref, n);

   // test output and constructor correct assignment of values
   cout << "Testing constructor and output " << endl;
   cout << "exp Taylor series : ";
   f_exp.print();
   cout << "cos Taylor series : ";
   f_cos.print();
   cout << "geom Taylor series : ";
   f_geom.print();
   cout << "Poly Taylor series : ";
   f_poly.print();

   // test () operator
   double x;
   cout << "Testing () operator " << endl;
   cout << "First try, input a value for x1  " << endl;
   cin >> x;
   cout << "f_exp(x1) : " << f_exp(x) << endl;
   cout << "numerical value and ratio for exp(x1)  " << exp(x) << " " << f_exp(x) / exp(x) << endl;
   cout << "f_cos(x1) : " << f_cos(x) << endl;
   cout << "numerical value and ratio for cos(x1)  " << cos(x) << " " << f_cos(x) / cos(x) << endl;
   cout << "f_geom(x1) : " << f_geom(x) << endl;
   cout << "numerical value and ratio for 1/(1-x1)  " << 1. / (1. - x) << " " << f_geom(x) * (1 - x) << endl;
   cout << "f_poly(x1) : " << f_poly(x) << endl;
   cout << "numerical value and ratio for poly(x1)  " << -1. + 6. * x - 2. * x * x * x + 9. * pow(x, 4) + 3. * pow(x, 5) << " "
        << f_poly(x) / (-1. + 6. * x - 2. * x * x * x + 9. * pow(x, 4) + 3. * pow(x, 5)) << endl;

   cout << "Second try, input a value for x2  " << endl;
   cin >> x;
   cout << "f_exp(x2) : " << f_exp(x) << endl;
   cout << "numerical value and ratio for exp(x2)  " << exp(x) << " " << f_exp(x) / exp(x) << endl;
   cout << "f_cos(x2) : " << f_cos(x) << endl;
   cout << "numerical value and ratio for cos(x2)  " << cos(x) << " " << f_cos(x) / cos(x) << endl;
   cout << "f_geom(x2) : " << f_geom(x) << endl;
   cout << "numerical value and ratio for 1/(1-x2)  " << 1. / (1. - x) << " " << f_geom(x) * (1 - x) << endl;
   cout << "f_poly(x2) : " << f_poly(x) << endl;
   cout << "numerical value and ratio for poly(x2)  " << -1. + 6. * x - 2. * x * x * x + 9. * pow(x, 4) + 3. * pow(x, 5) << " "
        << f_poly(x) / (-1. + 6. * x - 2. * x * x * x + 9. * pow(x, 4) + 3. * pow(x, 5)) << endl;

   // test addition operator. Use of default copy constructor to avoid use of forbidden assignment operator with
   // already existing instance of FuncTaylor
   cout << "Testing addition " << endl;
   FuncTaylor f_cos_p_exp(f_cos + f_exp_ref);
   cout << "f_cos+f_exp : ";
   f_cos_p_exp.print();
   cout << "f_cos_p_exp(x2) : " << f_cos_p_exp(x) << endl;
   cout << "numerical value and ratio for cos(x2)+exp(x2)  " << cos(x) + exp(x) << " " << f_cos_p_exp(x) / (cos(x) + exp(x))
        << endl;

   // test subtraction
   cout << "Testing subtraction " << endl;
   FuncTaylor f_poly_m_exp(f_poly_ref - f_exp_ref);
   cout << "f_poly-f_exp : ";
   f_poly_m_exp.print();
   cout << "f_poly_m_exp(x2) : " << f_poly_m_exp(x) << endl;
   cout << "numerical value and ratio for -1+6x2-2x2^3+9x2^4+3x2^5-exp(x2)  "
        << -1. + 6. * x - 2. * x * x * x + 9. * pow(x, 4) + 3. * pow(x, 5) - exp(x) << " "
        << f_poly_m_exp(x) / (-1. + 6. * x - 2. * x * x * x + 9. * pow(x, 4) + 3. * pow(x, 5) - exp(x)) << endl;

   // test derivation
   FuncTaylor f_sin(f_cos.deriv(3u));
   cout << "f_sin from derivation : ";
   f_sin.print();
   cout << "f_sin(x2) : " << f_sin(x) << endl;
   cout << "numerical value and ratio for sin(x2)  " << sin(x) << " " << f_sin(x) / sin(x) << endl;

   FuncTaylor f_null(f_poly.deriv(n + 2));
   cout << "f_null from derivation : ";
   f_null.print();
   cout << "f_null(x2) : " << f_null(x) << endl;

   FuncTaylor f_cst(f_poly.deriv(n));
   cout << "f_cst from derivation : ";
   f_cst.print();
   cout << "f_cst(x2) : " << f_cst(x) << endl;

   // test product
   FuncTaylor f_prod_cos_exp(f_cos * f_exp_ref);
   cout << "f_cos*f_exp : ";
   f_prod_cos_exp.print();
   cout << "f_prod_cos_exp(x2) : " << f_prod_cos_exp(x) << endl;
   cout << "numerical value and ratio for cos(x2)*exp(x2)  " << cos(x) * exp(x) << " " << f_prod_cos_exp(x) / (cos(x) * exp(x))
        << endl;

   FuncTaylor f_prod_geom_exp(f_geom_ref * f_exp_ref);
   cout << "f_geomxf_exp : ";
   f_prod_geom_exp.print();
   cout << "f_prod_geom_exp(x2) : " << f_prod_geom_exp(x) << endl;
   cout << "numerical value and ratio for exp(x2)/(1-x2)  " << exp(x) / (1 - x) << " " << f_prod_geom_exp(x) * (1 - x) / exp(x)
        << endl;

   // test "complex"  operation
   cout << "Testing a complex operation " << endl;
   FuncTaylor f_tc(f_geom_ref - f_cos * f_exp_ref);
   cout << "f_tc : ";
   f_tc.print();
   cout << "f_tc(x2) : " << f_tc(x) << endl;
   cout << "numerical value and ratio for f_tc(x2)  " << 1 / (1 - x) - cos(x) * exp(x) << " "
        << f_tc(x) / (1 / (1 - x) - cos(x) * exp(x)) << endl;

   // test getPoint method
   cout << "point_tc : " << f_tc.getPoint() << endl;

   // test exception for problematic cases
   cout << "Testing exceptions from problematic cases " << endl ;
   bool catched=false;
   try
   {
       FuncTaylor pb(f_cos*f_exp);
       FuncTaylor pbp(f_cos+f_geom);
       FuncTaylor pbm(f_cos-f_poly);
   }
   catch (int e)
   {
        cout << "Error1 : catching exception for note unique a in operation"<<endl;
        catched=true;
   }
   if (!catched)
        cout << "Error1 : not catching exception for note unique a in operation"<<endl;

   try
   {
       f_cos=f_exp;
   }
   catch (int e)
   {
        cout << "Error2 : catching exception for = usage"<<endl;
   }
   return 0;
}
