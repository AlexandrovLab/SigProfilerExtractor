CXX      := g++
NVCXX    := nvcc
CXXSTD   := -std=c++11
CXXFLAGS := $(CXXSTD) -O3 -shared -arch=sm_70 # This last flag makes the code exclusive for NVIDIA V100
NVOPTS   := --compiler-options -fPIC 
INCLUDES := $(shell python -m pybind11 --includes)
PYSUFFIX := $(shell python3-config --extension-suffix)

all: cuda_nmf_single$(PYSUFFIX) cuda_nmf_double$(PYSUFFIX)

%$(PYSUFFIX): %.cu
	$(NVCXX) $(CXXFLAGS) $(NVOPTS) $(INCLUDES) $< -o $@ -lnvToolsExt -lcublas