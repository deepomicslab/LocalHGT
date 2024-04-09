CC ?= g++
CFLAGS = $(CXXFLAGS) -std=c++11 -pthread 
LDFLAGS =

extract_ref: ./src/extract_ref_normal_peak.cpp
	$(CC) $(CFLAGS) ${LDFLAGS} -o scripts/extract_ref ./src/extract_ref_normal_peak.cpp
