include ../config/make.def

# Note that COMPILE is also defined in make.common and should
# be the same. We can't include make.common because it has a lot
# of other garbage. 
FCOMPILE = $(F77) -c $(F_INC) $(FFLAGS)

all: setparams 

# setparams creates an npbparam.h file for each benchmark 
# configuration. npbparams.h also contains info about how a benchmark
# was compiled and linked

setparams: setparams.cpp ../config/make.def
	$(UCC)  -L/home/titansarus/Documents/CMSC714/project/my/tbb-2020.3-lin/tbb/lib/intel64/gcc4.8 -Wl,-rpath,/home/titansarus/Documents/CMSC714/project/my/tbb-2020.3-lin/tbb/lib/intel64/gcc4.8 -I/home/titansarus/Documents/CMSC714/project/my/tbb-2020.3-lin/tbb/include/ -Wno-changes-meaning -ltbb -o setparams setparams.cpp 

clean: 
	-rm -f setparams setparams.hpp npbparams.hpp
	-rm -f *~ *.o

