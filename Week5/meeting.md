
[1] FFTW computes an unnormalized DFT. Thus, computing a forward followed by a backward transform (or vice versa) results in the original array scaled by n. For the definition of the DFT, see What FFTW Really Computes. 

* [https://www.fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html#Complex-One_002dDimensional-DFTs]

[2] Filtering algorithm and Random phase algorithm mainly differs from how to generate white noise


[3] simplification in week3 only validates for elastic normal contact, and may work for similar time-scale viscoelastic


[4] If we want to apply CMakeLists.txt in one project, we need first:

```bash
mkdir build
cd build/
cmake ..
make
```
**cmake** will generate *Makefile*, **make** will generate excutable file(*our_test*):

```bash
(base) zli@dionysos:~/Master-internship-repo-at-Institut-dAlembert/build$ cmake ..
-- Configuring done
-- Generating done
-- Build files have been written to: /gagarine/temporaires/zli/Master-internship-repo-at-Institut-dAlembert/build
(base) zli@dionysos:~/Master-internship-repo-at-Institut-dAlembert/build$ make
[ 50%] Building CXX object CMakeFiles/our_test.dir/hello_world_for_cpp_configure/test_for_fftw.cc.o
[100%] Linking CXX executable our_test
[100%] Built target our_test
(base) zli@dionysos:~/Master-internship-repo-at-Institut-dAlembert/build$ ls
CMakeCache.txt  CMakeFiles  cmake_install.cmake  compile_commands.json  Makefile  our_test
```

**ldd our_test** will check the state of linkage:

```bash
(base) zli@dionysos:~/Master-internship-repo-at-Institut-dAlembert/build$ ldd our_test 
        linux-vdso.so.1 (0x00007ffe764aa000)
        libfftw3.so.3 => /gagarine/temporaires/zli/fftw/lib/libfftw3.so.3 (0x00007f5b0ccd9000)
        libstdc++.so.6 => /lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007f5b0ca91000)
        libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007f5b0c869000)
        libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007f5b0c782000)
        /lib64/ld-linux-x86-64.so.2 (0x00007f5b0ce0a000)
        libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007f5b0c762000)
```

In *CMakeLists.txt*, there are some parameters that we can play with, *Debug* can collaborate with *gdb* for debugging, however, *Release* will give faster compiling speed.

```bash
//No help, variable specified on the command line.
CMAKE_BUILD_TYPE:STRING=Debug
```
In terminal, we can also use the following command to *explicitly* define:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake -DCMAKE_BUILD_TYPE=Debug ..
```
where *-D* for 'define'

[5] try to define a class for accessing 2-D array in C++(guard block for head file)
    compare/call library mdspan: [https://github.com/kokkos/mdspan](require C++ 23)
                                 [https://godbolt.org/]

    we just get *h*, try to replace the *-r**2/(2*R)* one for real contact, better with C++ solver.

[6] understand the rest script of Yastrebov and try to implement viscoelastic part