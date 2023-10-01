#RunMorphologies : $(SRC)
#icc -I $(BOOST_ROOT)/include $(SRC) -o ../2_run/RunMorphologies -openmp -O2  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -std=c++11
#
#	clang++ -I $(BOOST_ROOT) $(SRC) -o ../2_run/RunMorphologies -std=c++11 -O2 -lpthread  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -fopenmp -liomp5 
BOOST_ROOT=/usr/local/include/boost

SRC = morphologies.cpp shake.cpp initialization.cpp

RunMorphologies : $(SRC)
	icc -I $(BOOST_ROOT) $(SRC) -o ../2_run/RunMorphologies -std=c++11 -O2