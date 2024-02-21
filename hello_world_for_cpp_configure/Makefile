# compilation flag for optimized (Release mode)
CXXFLAGS = -O3 -fopenmp -DNDEBUG -std=c++17 -Wfatal-errors 
# compilation flag for debug mode. To use comment the above line and uncomment the line below
#CXXFLAGS = -g -DNDEBUG -std=c++11 -Wfatal-errors

# BLASLIB is the flag to link with your BLAS implementation. below is the default: the system bls lib. 
BLASLIB = -lblas

# FFTW is the flag to link with your FFTW implementation. 
FFTW = -lfftw

# To use an other blas complete and uncomment the line below, set for openblas
# You can get the openblas source from here : 
# www.openblas.net, get the tar.gz
#  copy it in your directory 
#  extract it : tar -zxvf OpenBLAS-0.3.26.tar.gz
#  compile it : cd OpenBLAS-0.3.26
#               make
#  When it is done you should have the file libopenblas.a, the openblas library

OPENBLASLIB = /gagarine/temporaires/zli/OpenBLAS-0.3.26/libopenblas.a -lpthread 

FFTW_DIR=/gagarine/temporaires/zli/fftw

# this list all the executable that will be build when typing make
#all: dmatrix_dense_examples test_prod_blas test_prod_openblas
all:


# this are the rule to compile the library file that Igave you (mmio.cpp and dmatrix_denseCM.cpp)
dmatrix_denseCM.o : Makefile dmatrix_denseCM.cpp mmio.o mmio.h dmatrix_denseCM.h
	$(CXX) $(CXXFLAGS) -c dmatrix_denseCM.cpp

mmio.o :  Makefile mmio.cpp mmio.h
	$(CXX) $(CXXFLAGS) -c mmio.cpp	

dmatrix_dense_examples.o : Makefile dmatrix_dense_examples.cpp dmatrix_denseCM.h mmio.h
	$(CXX) $(CXXFLAGS) -c dmatrix_dense_examples.cpp

dmatrix_dense_examples : Makefile dmatrix_dense_examples.o dmatrix_denseCM.o mmio.o
	$(CXX) $(CXXFLAGS) dmatrix_dense_examples.o dmatrix_denseCM.o mmio.o $(OPENBLASLIB) -o dmatrix_dense_examples


main.o : Makefile main.cpp dmatrix_denseCM.h mmio.h
	$(CXX) $(CXXFLAGS) -c main.cpp


test_prod_blas : Makefile main.o dmatrix_denseCM.o mmio.o
	$(CXX) $(CXXFLAGS) main.o dmatrix_denseCM.o mmio.o $(BLASLIB) -o test_prod_blas

test_prod_openblas : Makefile main.o dmatrix_denseCM.o mmio.o
	$(CXX) $(CXXFLAGS) main.o dmatrix_denseCM.o mmio.o $(OPENBLASLIB) -o test_prod_openblas

clean:
	rm *.o
	rm dmatrix_dense_examples test_prod_blas test_prod_openblas


