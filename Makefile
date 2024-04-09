CXX := g++
CXXFLAGS := $(CXXFLAGS) -std=c++11 -pthread 
extract_ref: ./src/extract_ref_normal_peak.cpp
	$(CXX) $(CXXFLAGS) -o scripts/extract_ref ./src/extract_ref_normal_peak.cpp
