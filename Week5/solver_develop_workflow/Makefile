# compilation flag for optimized (Release mode)
CXXFLAGS = -O3 -fopenmp -DNDEBUG -std=c++17 -Wfatal-errors 

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
all: random_roughness

random_roughness : Makefile test_for_fftw.cc
	$(CXX) $(CXXFLAGS) -o random_roughness test_for_fftw.cc -I$(FFTW_DIR)/include -L$(FFTW_DIR)/lib -lfftw3

clean:
	rm *.o
	rm random_roughness