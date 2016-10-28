BOOST_PATH=${HOME}/C++/boost_1_56_0

INCLUDES=-I${BOOST_PATH} 
LINKS=-L${BOOST_PATH}/lib 

LDFLAGS= -static -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lboost_serialization

CPP=g++
FF=ifort                                                                                                                                                    
FFLAGS=-O3 -m64 -openmp
CPPFLAGS= -O3 -m64 -openmp -DNDEBUG
#LDFLAGS=  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

#CPP=gcc
#CPPFLAGS= -O3 -m64 -fopenmp -DNDEBUG -fPIC

all:lapack_routines.o blasroutines.o helperfunctions.o
%.o: %.cpp 
	@echo 'Building $<'
	$(CPP) $(CPPFLAGS) $(INCLUDES) $(LINKS) -c $< -o $@ $(LDFLAGS)
	@echo 'Finished building $@'
	@echo ' '
clean:
	rm -f  *.o
