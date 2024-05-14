Since we are contributing [Tammas](https://gitlab.com/tamaas/tamaas), we first need to figure out the structure of tamaas works.

### src

#### model

##### model class

Model containing pressure and displacement. This class is a container for the model fields. It is supposed to be dimension agnostic, hence the GridBase members.

`model.hh` contains getDiscretization() function, we can use this function for *n* and *m* variable of our python script.





#### solvers

##### polonsky_keer_rey










##### Advantages of Using `std::vector` for User Input
1. **Familiarity and Accessibility**: `std::vector` is a standard C++ container, widely recognized and used by programmers. It offers dynamic size adjustment, direct access to elements, and is compatible with a wide range of standard library functions, making it an excellent choice for user-facing interfaces.

2. **Flexibility**: Users can easily construct, modify, and pass `std::vector` objects around in their applications. It supports a variety of operations and can be resized dynamically, which is ideal for handling user input where the size or content of data may not be known in advance.

3. **Safety and Performance**: `std::vector` manages its own memory and ensures that data is contiguous in memory, which optimizes performance especially for large data sets. This is crucial when user inputs are directly used in computationally intensive tasks.



