I want to add C++ now, there could be more attemptions that could output excutable files with less use, is any gitignore file necessary?

p.s. I get used to use bash mode for compiling, but debuging with vs code could be interesting, a little bit configuration has started in .vscode file folder.

I learned C++ 98, any advise for version? Let's have a look at .vscode file folder.

### about FFTW
I tried to install FFTW library, the constraint is that I don't have the access to install it in *root*. I install it in 

```bash
/gagarine/temporaires/zli/fftw
```
 with instructions[2], and tried to call fftw library in 

```bash
/gagarine/temporaires/zli/Master-internship-repo-at-Institut-dAlembert.
```

In *.bashrc*, I define

```bash
export FFTW_DIR=/gagarine/temporaires/zli/fftw
```

so to compile the cpp file, I have to use a long command

```bash
g++ -o fftw test_for_fftw.cc -I$FFTW_DIR/include -L$FFTW_DIR/lib -lfftw3
```
**The command can work, and gived me reference output.** But is there any way to simplify my command in terminal? with alias command? reinstall?

Besides, I didn't try mpi/openmp installation for FFTW:

```bash
./configure --prefix=/home/xxx/usr/fftw --enable-mpi --enable-openmp --enable-threads --enable-shared MPICC=mpicc CC=gcc F77=gfortran
```

Here is complete command for my installation

```bash
tar -xzvf fftw-3.3.10.tar.gz 
cd fftw-3.3.10/
./configure --prefix=/gagarine/temporaires/zli/fftw
make
make install
```
The lab computer system seems to have OpenMP:

```bash
(base) zli@dionysos:~$ gcc --version
gcc (Ubuntu 11.3.0-1ubuntu1~22.04) 11.3.0
Copyright (C) 2021 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
(base) zli@dionysos:~$ echo |cpp -fopenmp -dM |grep -i open
#define _OPENMP 201511
```

But MPI seems not to be configured successfully, it is different form my own Ububtu 22.04.

```bash
Configuration error:
Votre fichier de configuration comporte une erreur de programmation : 

Traceback (most recent call last):
  File "/gagarine/temporaires/zli/anaconda3/lib/python3.11/site-packages/sphinx/config.py", line 343, in eval_config_file
    exec(code, namespace)
  File "/gagarine/temporaires/zli/openmpi-5.0.2/docs/conf.py", line 166, in <module>
    import sphinx_rtd_theme
ModuleNotFoundError: No module named 'sphinx_rtd_theme'

make[1]: *** [Makefile:2712 : _build/man/ompi-wrapper-compiler.1] Erreur 2
make[1] : on quitte le répertoire « /gagarine/temporaires/zli/openmpi-5.0.2/docs »
make: *** [Makefile:1533 : all-recursive] Erreur 1
```

## Update 22/02/2024

I use CMake extension on the left bar, and it automatic generate a build folder that contains a series of CMake files(look powerful).

I also tried Makefile that I am familier with in ECN couse DDIS. Luckily, Makefile always work, it seems more suitable for large project, especially for linking excutable files in different steps.

So, are Makefile, CMake, CMakelists.txt different solutions, which one is better? I want to handle them all.

### Reference:

[1] https://code.visualstudio.com/docs/cpp/config-linux

[2] https://zhuanlan.zhihu.com/p/600161033

[3] 