#include "Vector2.h"
#include <cstring>
#include <new>
#include <cmath>

using namespace std;

Vector::Vector(const int & size_)
    : size(size_),values(NULL)
{

    // Here observe that if size<1 the second test is not done. It preserve trying allocating with
    // a 0 or negative value.
    if (size < 1 || !( values = new (nothrow) double[size] ))
    {
        cout<<"Error in Vector constructor : size is inappropriate"<<endl;
        if (size < 1) cout<<"Size parameter should be greater than 0"<<endl;
        else cout<<"Size is to big. Memory exhausted"<<endl;
        throw -size;
    }

}
Vector::Vector(const Vector & input)
    : size(input.size),values(NULL)
{
    if ( !( values = new (nothrow) double[size] ) )
    {
        cout<<"In copy constructor memory is exhausted"<<endl;
        throw -size;
    }
    // memcpy : This is a way to copy information stored at a specific address in memory to an other address.
    // This information can be a double a int or what so ever. But it also can be a set of double like here.
    // See man page for further information
    memcpy((void *) values,(void *) input.values,sizeof( double )*size);
}
Vector::~Vector()
{
    // why a test as constructor are checking things? Because we throw a exception when there is problem in constructor
    // so we don't know if this exception will not be cached. If yes program may continue with unallocated memory for this
    // instance. And as soon as the program jump out of the context where this instance have been created it will call
    // this method and delete should not be done as nothing was allocated. By testing if values is not null we prevent
    // bug.
    // This also explain why in constructor the first thing which is done is to set values to NULL.
    if (values) delete [] values;
    values = NULL;
}
Vector & Vector::operator=(const Vector & v)
{
    const int min_size = ( size > v.size ) ? v.size : size;
    memcpy((void *) values,(void *) v.values,sizeof( double )*min_size);
    // for debuging purpose only
    if (size != v.size)
        cout<<"warning: = operator on vector of different size ! Smallest part treated"<<endl;
    return *this;
}
double Vector::operator()(const int & i) const
{
    if ( i > -1 && i < size )
        return values[i];
    else
    {
        cout<<"error: incorrect index ! "<<i<<endl;
        throw -3;
    }
}
double & Vector::operator()(const int & i)
{
    if ( i > -1 && i < size )
        return values[i];
    else
    {
        cout<<"error: incorrect index ! "<<i<<endl;
        throw -4;
    }
}
Vector & Vector::operator*=(const double & a)
{
    for (int i = 0; i < size; ++i) values[i] *= a;
    return *this;
}
Vector & Vector::operator+=(const Vector & v)
{
    const int min_size = ( size > v.size ) ? v.size : size;
    for (int i = 0; i < min_size; ++i) values[i] += v.values[i];
    return *this;
}
Vector & Vector::operator-=(const Vector & v)
{
    const int min_size = ( size > v.size ) ? v.size : size;
    for (int i = 0; i < min_size; ++i) values[i] -= v.values[i];
    return *this;
}
Vector Vector::operator*(const double & a) const
{
    // without use of *= operator ; see below for alternative
    Vector v(size);
    for (int i = 0; i < size; ++i) v.values[i] = values[i]*a;
    return v;
}
Vector Vector::operator+(const Vector & v) const
{
    // this implementation use the += operator. From a performance point of view some may
    // say that using the copy constructor add the work of copying values of left hand side terms in
    // newly created memory. If we had put a full implementation here it would have been direct. In
    // a new vector we set by a loop directly the resulting sum.
    // Well it's true but the memcpy is rather efficient and the data are certainly mostly in the cache afterward.
    //
    // But the interesting aspect of this design is that you implement only once the for loop (in the += operator)
    // and this is the important fact. Only one implementation of the same thing avoid buggs.
    // Imagine that operator where more complicated, let say 20 lines of codes. Use of this technique (when possible)
    // avoid writing 40 lines of codes. And if bug have to be corrected you only need to correct one piece of code not
    // two.
    //
    Vector res_vect(*this);
    res_vect += v;
    return res_vect;
}
Vector Vector::operator-(const Vector & v) const
{
    Vector res_vect(*this);
    res_vect -= v;
    return res_vect;
}
double Vector::operator*(const Vector & v) const
{
    double res = 0.;
    if ( v.size != size )
    {
        cout<<"error:  for * operator vector must be of the same size !"<<endl;
        throw -5;
    }
    else
    {
        for (int i = 0; i < size; ++i) res += values[i]*v.values[i];
    }
    return res;
}
double Vector::norm(void) const
{
    // Again here we use * operator implemented above to avoid duplication
    // if we decide to us blas or manual unrolled loop to improve
    // performance of * operator, this norm operator will have also the
    // benefit of any improvement.
    double res = ( *this )*( *this );
    return sqrt(res);
}

int Vector::getSize(void) const
{
    return size;
}
std::ostream & operator << (std::ostream & ofs, const Vector & v)
{
    int i,size = v.size;
    if (size)
    {
        ofs << "[";
        ofs << v.values[0];
        for (i = 1; i < size; ++i)
        {
            ofs<<", "<< v.values[i];
        }
        ofs << "]";
    }
    else
        ofs <<endl<<" Empty vector ";
    return ofs;
}
