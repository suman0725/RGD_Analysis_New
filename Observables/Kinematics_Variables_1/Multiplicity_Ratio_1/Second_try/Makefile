CC=g++
CFLAGS=-std=c++17

# Include paths
INCLUDES=-I/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/include/hipo4 \
         -I/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.6/include \
         $(shell root-config --cflags)

# Linker flags
LDFLAGS=-L/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/lib \
         -L/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.6/lib \
         $(shell root-config --libs)

# Libraries to link against (NOTE the -lEG added!)
#TDatabasePDG, TParticlePDG, and related classes are defined in ROOT’s libEG.so (EG stands for “Event Generator”).

#root-config --libs does not include -lEG by default, which is why you had to add it manually.
LIBS=-lhipo4 -lccdb -lEG -pthread -lm -ldl

# Rule to compile any .cpp file to an executable with the same name
%: %.cpp
	$(CC) $(CFLAGS) -o $@ $< $(INCLUDES) $(LDFLAGS) $(LIBS)

# Clean rule
clean:
	rm -f *.o *~ core
