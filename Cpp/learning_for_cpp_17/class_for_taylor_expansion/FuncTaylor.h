#ifndef FuncTaylor_lab_H
#define FuncTaylor_lab_H
#include "Vector2.h"

// The FuncTaylor class stores a function as a Taylor expansion computed at a given point 'a' for a given order.
// The class use an implicit storage: Only the derivative of the function at point 'a' is stored by the class, considering
// that the 1/i! part can be computed on the fly when needed. 
// To evaluated the function represented by this class at a given point 'x', the point 'a' must also be stored.
// All computation are performed in real arithmetic (represented by doubles).
// Only the basic arithmetic around  this function concept is done here. 
// This is mostly an exercise to manipulate  member variable of class type. 
// 

class FuncTaylor {
    public:
     // Constructor member
     // Instantiation of a FuncTaylor object to represent an 'F' function.
     // input :
     //     fa_ = derivatives of the function 'F' at point 'a' . They are stored in ascending order of derivation in a continuous
     //     memory array pointed by this argument. 
     //     a_ = point at which 'F' and its derivatives are computed.
     //     order = order at which the expansion is computed.
     //
     //  There must be  order+1 values available in array pointed by fa_.
     FuncTaylor(const double* fa_, double a_, int order);

     // access member
     // Prints the Taylor expansion terms stored by the class (only the independent part of (x-a)^i/i!), to the standard output.
     void print(void) const;
     // obtain the point at which Taylor expansion terms were computed based on user input at construction time
     double getPoint(void) const;

     // computing member
     // Evaluating the function at a given value.
     // input :
     //    value x at which function is evaluated
     // output :
     //    function value
     double operator()(const double& x) const;

     // operator member
     // Addition of 2 functions
     // input : instance invoking this method (left hand side of the operator) called below q
     //         instance given as argument (right hand side of the operator) called below p
     // output : a new instance of the class FuncTaylor, containing q+p
     // operator works only in between Taylor expansion computed at the same point
     FuncTaylor operator+(const FuncTaylor& p) const;

     // operator member
     // Subtraction of 2 functions
     // input : instance invoking this method (left hand side of the operator) called below q
     //         instance given as argument (right hand side of the operator) called below p
     // output : a new instance of the class FuncTaylor, containing q-p
     // operator works only in between Taylor expansion computed at the same point
     FuncTaylor operator-(const FuncTaylor& p) const;

     // operator member
     // Multiplication of 2 function
     // input : instance invoking this method (left hand side of the operator) called below q
     //         instance given as argument (right hand side of the operator) called below p
     // output : a new instance of the class FuncTaylor, containing q*p
     // operator works only in between Taylor expansion computed at the same point
     FuncTaylor operator*(const FuncTaylor& p) const;

     // computing member
     // return as a new instance the Nth derivative of the function
     FuncTaylor deriv(const unsigned int& nth) const;

     // for safety
     // To avoid use of this operator, implementation of a fake = operator throwing a error
     // Limitation of Vector class don't permit full implementation of this operator. As
     // mentioned it's a exercise
     FuncTaylor& operator=(const FuncTaylor& p);

    private:
        // data member to store 'a' used when evaluating function for a specific value
        double a;
        // data member of type Vector to store the value that multiply terms (x-a)^i/i! i=0,1,...,order.
        // At construction times it corresponds to F(a),F'(a),F''(a), .....
        Vector terms;
};

#endif
