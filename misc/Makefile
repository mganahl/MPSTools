BOOST_PATH=${HOME}/C++/boost_1_54_0

INCLUDES=-I${BOOST_PATH} 
LINKS= -L${BOOST_PATH}/lib 

LDFLAGS= -static -lBoostTools -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lboost_serialization
#LDFLAGS= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lboost_serialization

CPP=g++
FF=ifort                                                                                                                                                    
FFLAGS=-O3 -m64 -openmp
CPPFLAGS= -O3 -m64 -openmp -DNDEBUG -fPIC

#CPP=gcc
#CPPFLAGS= -O3 -m64 -fopenmp -DNDEBUG -fPIC

all:Utilities.o File.o
%.o: %.cpp 
	@echo 'Building $<'
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@ $(LDFLAGS)
	@echo 'Finished building $@'
	@echo ' '
clean:
	rm -f  *.o
