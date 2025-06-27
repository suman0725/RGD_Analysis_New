CC=g++
CFLAGS=-std=c++17  # Use C++17 instead of C++11
# Include paths
INCLUDES=-I/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/include/hipo4 \
         -I/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.6/include \
         $(shell root-config --cflags)  # Include ROOT headers

# Linker flags
LDFLAGS=-L/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/lib \
         -L/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.6/lib \
         $(shell root-config --libs)  # Link ROOT libraries

# Libraries to link against
LIBS=-lhipo4 -lccdb -pthread -lm -ldl

# Rule to compile any .cpp file to an executable with the same name
%: %.cpp
	$(CC) $(CFLAGS) -o $@ $< $(INCLUDES) $(LDFLAGS) $(LIBS)

# Clean up build files
clean:
	rm -f *.o *~ core
