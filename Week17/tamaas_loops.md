The critical part in our viscoelastic implementation for tamaas project is to deal with a for loop in Python, but actually in C++ solver tamaas, the 'loop' is user-defined, we only implement `solve` function as one time-step. 

```cpp
Loop::loop([alpha, beta, gamma](Real& h_new, const Real& surface, const Real& u, const  Real & maxwell) {
    h_new = alpha * surface - beta * u + maxwell;
  }, H_new, this->surface, U, M_maxwell);
```


```cpp
auto func = [](...) { };

double square(double x) {
  return x*x;
}


double alpha = 1.;

auto square = [alpha](double x) {
  return x*x + alpha;
};

square(3) == 10
```