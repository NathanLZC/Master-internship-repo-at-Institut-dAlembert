#ifndef LAB_REAL_H
#define LAB_REAL_H


class Real
{
    public:
        // Constructor with two arguments representing respectively the mantissa and exponent of
        // the real number we want to store in instance created.
        // A default zero value is used for both argument.
        // Use of such defaulted argument make this constructor usable with two, one or zero argument.
        // In this last case  it is equivalent in usage to a default constructor. But it force
        // members to be nullify which is not the case of a default constructor.
        Real(double m = 0., int e = 0);

        // display method
        void display() const;

        // addition operator
        Real operator+(const Real & rhs) const;

        // substraction operator
        Real operator-(const Real & rhs) const;

        // multiplication operator
        Real operator*(const Real & rhs) const;

        // division operator
        Real operator/(const Real & rhs) const;

        // power of n
        Real power(int n) const;

        // exponential
        Real exp() const;

    private:
        double m;
        int e;
};

#endif
