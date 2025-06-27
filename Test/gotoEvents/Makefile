CC=g++
CFLAGS=-std=c++11
INCLUDES=-I/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/include/hipo4
LDFLAGS=-L/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/lib
LIBS=-lhipo4 -llz4 $(shell root-config --cflags --glibs)

# Rule to compile any .cpp file to an executable with the same name
%: %.cpp
	$(CC) $(CFLAGS) -o $@ $< $(INCLUDES) $(LDFLAGS) $(LIBS)

# Clean up build files
clean:
	rm -f *.o *~ core
