#ifndef Vector_lab_H
#define Vector_lab_H
#include <iostream>

// Class Vector is defining a new type of object that deal with multidimensional Vector.
// Each instance of this class have it's own dimension and set of double values.
//
class Vector
{
    public:
        // Constructor member
        // Instantiate a Vector object using size argument to define the dimension of the vector
        // (i.e. the number of values this object will treat)
        // The object generated don't have values set after construction. They are undefined.
        // Incorrect size will make this method throwing an exception.
        // By incorrect we mean size<1 or size so big that you exhaust computer memory.
        Vector(const int & size_);

        // Constructor member
        // Instantiate a Vector object by using a other instance (copy constructor)
        // New instance will have the same dimension and values as input parameter
        Vector(const Vector &input);

        // Destructor
        ~Vector();

        // operator member
        // equal operator
        // input : instance invoking this method (left hand side of the operator) called below u
        //         instance given as argument (right hand side of the operator) called below v
        // output : a reference to u containing for the first common index the value of v.
        // By "first common index" we are talking of index starting at 0 up to the minimum of dimension v and
        // dimension u.
        Vector & operator=(const Vector & v);

        // operator member
        // give the ith value of a instance without changing this instance.
        // It's invocation is  in a const context (const instance)
        // input : instance invoking this method (left hand side of the operator) called below u
        //         query index i given in parentheses
        // output : double corresponding to index i of u
        double operator()(const int &i) const;

        // operator member
        // give the reference to the ith value of a instance
        // input : instance invoking this method (left hand side of the operator) called below u
        //         query index i given in parentheses
        // output : reference to double corresponding to index i of u
        double &operator()(const int &i);

        // operator member
        // scaling of a vector
        // input : instance invoking this method (left hand side of the operator) called below u
        //         scaling value a
        // output : a reference to u containing the product of u(i) by a, i in {0,...size-1}.
        Vector &operator*=(const double & a);

        // operator member
        // addition of a vector to an other
        // input : instance invoking this method (left hand side of the operator) called below u
        //         instance given as argument (right hand side of the operator) called below v
        // output : a reference to u containing for the first common index the sum of u and v.
        // By "first common index" we are talking of index starting at 0 up to the minimum of dimension v and
        // dimension u.
        Vector &operator+=(const Vector & v);

        // operator member
        // subtraction of a vector to an other
        // same as += except that arithmetic operation is subtraction.
        Vector &operator-=(const Vector & v);

        // operator member
        // Result of the multiplication of a vector by value.
        // input : instance invoking this method (left hand side of the operator) called below u
        //         scaling value a
        // output : a new instance of the class Vector containing u scaled by a.
        Vector operator*(const double & a) const;

        // operator member
        // addition of 2 vectors.
        // input : instance invoking this method (left hand side of the operator) called below u
        //         instance given as argument (right hand side of the operator) called below v
        // output : a new instance of the class Vector,initiated as a copy of u, and
        // containing the sum of first common index of u and v.
        // By "first common index" we are talking of index starting at 0 up to the minimum of dimension v and
        // dimension u.
        Vector operator + (const Vector & v) const;

        // operator member
        // subtraction of 2 vectors.
        // input : instance invoking this method (left hand side of the operator) called below u
        //         instance given as argument (right hand side of the operator) called below v
        // output : a new instance of the class Vector,initiated as a copy of u, and
        // containing the difference of first common index of u and v.
        // By "first common index" we are talking of index starting at 0 up to the minimum of dimension v and
        // dimension u.
        Vector operator - (const Vector & v) const;

        // operator member
        // dot product of 2 vector.
        // input : instance invoking this method (left hand side of the operator) called below u
        //         instance given as argument (right hand side of the operator) called below v
        // output : a double corresponding to the dot product
        // if you declare
        // double x=u*v;
        // then
        // x=Sum(u(i)*v(i))
        //    i in {0,...size-1}
        // If u and v don't have the same size, a exception is throw as it is not considered to be possible here.
        // This operator is returning a double corresponding to the dot product
        double operator * (const Vector & v) const;

        // computing member
        // Euclidean norm of the vector.
        // No argument used as it apply on instance which call this method
        // This member function is returning a double corresponding to the norm of the vector
        double norm(void) const;

        // access member
        //  give the dimension of the instance of the vector
        int getSize(void) const;


        // friend function to the class
        // to simplify output overloading of << operator of the class ostream to output instance of Vector class.
        friend std::ostream & operator << (std::ostream &, const Vector &);

    private:
        int size;
        double *values;
};

#endif
