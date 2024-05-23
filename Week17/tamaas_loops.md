The critical part in our viscoelastic implementation for [tamaas](https://gitlab.com/tamaas/tamaas) project is to deal with a "for loop" in Python, but actually in C++ solver tamaas, the solving loops is user-defined(iteration times, tolerance), we only implement `solve` function as one time-step of "for loop". 


#### One line to update surface profile $H_{t+\Delta t}$ should be noticed:

$$
\mathbb{H}^\prime = (\underbrace{G_{\infty}+\tilde{G}}_\alpha) \mathbb{H}-\underbrace{\tilde{G}}_\beta \mathbb{U_t}+\sum_k \underbrace{\frac{\tau^k}{\tau^k+\Delta t}}_{\gamma^k} \mathbb{M}_t^k
$$

One way to deal with these three coefficients is to use **Lambda function**, here one possible implementation can be given as:

```cpp
Loop::loop([alpha, beta, gamma](Real& h_new, const Real& surface, const Real& u, const  Real & maxwell) {
    h_new = alpha * surface - beta * u + maxwell;
  }, H_new, this->surface, U, M_maxwell);
```

Generally speaking, we have a **Lambda function** in this format:

```cpp
auto func = [](...) { };
```

We can use the following scripts as an example:

- Define `square` function in different namespace:
  
```cpp
#include <iostream>

namespace global_scope {
    double square(double x) {
        return x * x;
    }
}

namespace lambda_scope {
    double alpha = 1.;

    auto square = [alpha](double x) {
        return x * x + alpha;
    };
}

int main()
{
    std::cout << global_scope::square(3) << std::endl;  
    std::cout << lambda_scope::square(1) << std::endl;  

    return 0;
}
```

- We can also use reference since we plan to use **Lambda function** capture resuable parameters. For large objects or variables that need to be modified, it may be more appropriate to use reference capture.

```cpp
namespace lambda_scope {
    double alpha = 1.;

    auto square = [&alpha](double x) {  
        return x * x + alpha;
    };
}
```

### A good example for reference capture:

```cpp
int x = 4;
 
auto y = [&r = x, x = x + 1]() -> int
{
    r += 2;
    return x * x;
}(); // updates ::x to 6 and initializes y to 25.
```

This C++ code involves a lambda expression with reference and value captures, which results in updating an external variable and returning a calculated value.

First, the code initializes an integer variable `x`:

```cpp
int x = 4;
```

Then, a lambda expression is defined and immediately invoked, with its result assigned to `y`:

```cpp
auto y = [&r = x, x = x + 1]() -> int
{
    r += 2;
    return x * x;
}();
```

#### Capture List
The capture list `&r = x, x = x + 1` does the following:
- `&r = x`: Captures the external variable `x` by reference and names it `r`.
- `x = x + 1`: Captures the external variable `x` by value, initializing the captured value to `x + 1`, which is `5`.

Therefore, within the lambda:
- `r` is a reference to the external `x`.
- `x` is a local variable with a value of `5`.

#### Lambda Body
The body of the lambda expression is:

```cpp
{
    r += 2;
    return x * x;
}
```

- `r += 2`: Since `r` is a reference to the external `x`, this updates the external `x` to `4 + 2 = 6`.
- `return x * x`: This returns the square of the local variable `x`, which is `5 * 5 = 25`.

#### Result
The lambda is immediately invoked, and its return value is assigned to `y`.

Putting it all together:
1. The external variable `x` is updated to `6` (via the reference capture `r`).
2. The lambda returns `25` (the square of the local variable `x`), which is assigned to `y`.

#### Summary
- The external variable `x` is updated to `6`.
- `y` is initialized to `25`.
