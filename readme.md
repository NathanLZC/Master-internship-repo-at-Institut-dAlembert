# Master Thesis Project

## Fluage et surfaces rugueuses : évolution du contact entre matériaux viscoélastiques

We are aiming to design a solver for evolution of the contact between two random rough viscoelastic surfaces. If possible, we want to continue developing with wavelet method. This repo records our development workflow.

### Key word: FFT, Boundary element method, Viscoelastic, contact mechanics

## Getting Started

These instructions will help you get a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

This project is mainly based on Python and C++, the notebooks in */Weekn* contain theory background for contact mechanics. The python version is:

```bash
Python 3.10.6 64-bit
```

### Installing

1. **Clone the repo**

```bash
git clone https://github.com/NathanLZC/Master-internship-repo-at-Institut-dAlembert
```

2. **Navigate to the project directory**

```bash
cd Master-internship-repo-at-Institut-dAlembert
```

3. **Install the required packages**

Make sure you have `pip` installed and then run:

```bash
pip install -r requirements.txt
```

## Running the tests

### For .py and .ipynb files, we strongly recommand to use VS code.
### For .cc/.cpp files, we rely on Makefile and CMake


## Requirement libraries
### FFTW

### Eigen

With CMakeLists.txt, we can specify the location of the Eigen library. On the system that we have full access to, we can also use the following way to configure a global environment for Eigen:

First, download the required version of Eigen from Eigen's [official website](http://eigen.tuxfamily.org/) or its [GitLab repository](https://gitlab.com/libeigen/eigen). Take `eigen-3.4.0.tar.gz` as an example:

Then use the following code to decompress and create a symbolic link:
```bash
tar -xzf eigen-3.4.0.tar.gz -C /home/username
sudo ln -s /path/to/eigen-3.4.0/Eigen /usr/local/include/Eigen
```

### mdspan for C++23

## Streamlit


## Built With

* [Python](https://www.python.org/downloads/release/python-3106/) - The main language used
* [C++](https://cplusplus.com/) - The main language used
* [FFTW](https://www.fftw.org/) - A fast, free C FFT library; includes real-complex, multidimensional, and parallel transforms.
* [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) - A high-level C++ library of template headers for linear algebra, matrix and vector operations, geometrical transformations, numerical solvers and related algorithms.

## Authors

* **Zichen Li** - *Initial work* - [NathanLZC](https://github.com/YourUsername)

## License

This project is licensed under the CC0-1.0 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
