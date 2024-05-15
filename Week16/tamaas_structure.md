Since we are contributing [Tammas](https://gitlab.com/tamaas/tamaas), we first need to figure out the structure of tamaas works.

 We can see `namespace tamaas{}` a lot in our project, in C++ programming, the namespace feature is highly important as it provides a way to organize identifiers (such as variables, functions, and classes) into separate logical groups. This helps to avoid naming conflicts and enhances code readability. The primary functions and advantages of namespaces include:

1. **Avoiding Naming Conflicts**:
   In large projects, different parts might use the same identifier names. Without namespaces, these identical names can lead to conflicts and confusion. 

2. **Organizing Code**:
   Namespaces help logically organize code, making it easier to maintain and understand. 

3. **Providing Aliases**:
   Namespaces allow the creation of aliases to simplify the use of long namespace names, improving code readability and usability. 

4. **Allowing Nested Namespaces**:
   Namespaces can be nested, which adds more hierarchical organization and logical structure to the code. 

5. **Using Anonymous Namespaces for Internal Linkage**:
   Anonymous namespaces restrict the visibility of their contents to the file in which they are declared, preventing unnecessary external linkage.

6. **Using Namespaces**:
   Namespaces can be used in three ways: fully qualified names, `using` declarations, and `using` directives.


The source code for tamaas solver is under `/src`:

### src

#### 1. core

This folder holds low-level staff for tammas project, we can find `cuda`, `fftw`, `loops`, which provide interface such that we can call third-party library. And this folder also contains a series of scripts like:

##### grid_base.hh

The GridBase class template not only provides functionality similar to an `array` but also includes support for parallelization through MPI and CUDA. 



#### 2. model

This folder contains the implementation of Green's function, material law and energy minimization. 

##### model class

Model containing pressure, displacement and operators for elasticity. This class is a container for the model fields. It is supposed to be dimension agnostic, hence the GridBase members.

`model.hh` contains `getDiscretization()` function, we can use this function for *n* and *m* variable of our python script. It also contains `getShearModulus()` function, since we implement a generalized Maxwell model in our codes, the return value of this function can represent $G_{\infty}$, which is shear modulus of pure elastic branch.

##### westergaard class

FFT and fourrier coefficient are well implemented in this class.


#### 3. solvers
Our implementation for viscoelasticity with generalized Maxwell model should derive from PolonskyKeerRey class, the critical function is `Real solve(std::vector<Real> target)`, which is first derived from base class `ContactSolver` by class `PolonskyKeerRey`.

##### polonsky_keer_rey

??
The following two types are for pressure:
```bash
protected:
  type variable_type, constraint_type;
```

Additionally, `primal` represents pressure, which is the paramter to be optimized, and `dual` represents gap(gradient).


Loop

target







#### 4. surface









