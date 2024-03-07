#ifndef _ITERATIVES_SOLVERS_
#define _ITERATIVES_SOLVERS_

#include <iostream>
#include <vector>

template < class V_type >
struct iter_solver_return_type{
  bool converged; // The itrative method will set thiss value to true if it converged.
  int nbiter;     // number of iteration to reach convergence
  double rnorm_res; // relative residual norm for the solution ( norm(Ax -b)/ norm(b))
  V_type x; // the solution
  std::vector<double > iter_rnorm_res; // a vector containing the relative residual for each iteration. filled up only if the parameter monitor is set to true when called the iterative solver. 
};

/// A short cut to print the result of the iterative method
template <class T>
std::ostream & operator << ( std::ostream & out, const iter_solver_return_type<T> & results){
  if (results.converged){
    out << "Success in " << results.nbiter << " Iteration. " << "res norm " << results.rnorm_res; 
  }
  else{
    out << "Failure after " << results.nbiter << " iteration " << "res norm " << results.rnorm_res << std::endl;
  }
  return out;
}









#endif