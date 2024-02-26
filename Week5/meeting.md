
[1] FFTW computes an unnormalized DFT. Thus, computing a forward followed by a backward transform (or vice versa) results in the original array scaled by n. For the definition of the DFT, see What FFTW Really Computes. 

* [https://www.fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html#Complex-One_002dDimensional-DFTs]

[2] Filtering algorithm and Random phase algorithm mainly differs from how to generate


[3]


[4] CMakeLists.txt

//No help, variable specified on the command line.
CMAKE_BUILD_TYPE:STRING=Debug

[5] try to define a class for accessing 2-D array in C++(guard block for head file)
    compare/call library mdspan: [https://github.com/kokkos/mdspan](require C++ 23)
                                 [https://godbolt.org/]


[6] understand the rest script of Yastrebov and try to implement viscoelastic part