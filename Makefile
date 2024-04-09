# extract_ref: ./src/extract_ref_normal_peak.cpp
# 	g++ -pthread -std=c++11 -o scripts/extract_ref ./src/extract_ref_normal_peak.cpp

CXX := g++
CXXFLAGS := $(CXXFLAGS) -pthread -std=c++11
extract_ref: ./src/extract_ref_normal_peak.cpp
	$(CXX) $(CXXFLAGS) -o scripts/extract_ref ./src/extract_ref_normal_peak.cpp
