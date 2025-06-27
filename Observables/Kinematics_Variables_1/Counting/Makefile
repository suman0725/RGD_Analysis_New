CC=g++
CFLAGS=-std=c++17
# Include paths
INCLUDES=-I/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/include/hipo4

# Linker flags
LDFLAGS=-L/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/lib

# Libraries to link against
LIBS=-lhipo4 -pthread -lm -ldl

# Target to build count_files_events
count_files_events: count_files_events.cpp
	$(CC) $(CFLAGS) -o $@ $< $(INCLUDES) $(LDFLAGS) $(LIBS)

# Clean up build files and executable
clean:
	rm -f *.o *~ core count_files_events