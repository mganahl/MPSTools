INTEL=${HOME}/intel
BOOST_PATH=${HOME}/C++/boost_1_56_0
MPSTOOLS=${HOME}/workspace/MPSTools/MPSTools

#CPP=icpc


#FF=ifort                                                                                                                                             
#FFLAGS=-O3 -m64 -fopenmp
#CPPFLAGS=-m64 -openmp  -g3
CPP=g++
CPPFLAGS=-O3 -m64 -openmp -DNDEBUG -std=c++11 
OBJ = ../blas/lapack_routines.o ../blas/blasroutines.o ../blas/helperfunctions.o ../Parser/ParameterParser.o ../misc/File.o ../misc/Utilities.o


INCLUDES=-I${BOOST_PATH}  -I${MPSTOOLS}/sourcefiles/blas -I${MPSTOOLS}/sourcefiles/Parser -I${MPSTOOLS}/sourcefiles/src -I${MPSTOOLS}/sourcefiles/misc
#INCLUDES=-I${BOOST_PATH}  -I${MPSTOOLS}/sourcefiles

LINKS=-L${INTEL}/mkl/lib/intel64 -L${BOOST_PATH}/lib -L${MPSTOOLS}/sourcefiles/misc -L${MPSTOOLS}/sourcefiles/Parser -L${MPSTOOLS}/sourcefiles/blas
#LINKS=-L${INTEL}/mkl/lib/intel64 -L${BOOST_PATH}/lib -L${MPSTOOLS}/sourcefiles

LDFLAGS=  -static -Wl,--start-group -ldl -lirc -lc -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lboost_filesystem -lboost_system -lboost_serialization -Wl,--end-group -lpthread -lm 
#DYNAMIC LINKING
#LDFLAGS= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lboost_filesystem -lboost_system -lpthread -lboost_serialization

bin/%.exe: %.cpp
	@echo 'Building $<'
	$(CPP) $(CPPFLAGS) $(INCLUDES) $(LINKS) -o $@ $<  $(OBJ) $(LDFLAGS)
	@echo 'Finished building $@'
	@echo ' '

clean:
	rm bin/*.exe
