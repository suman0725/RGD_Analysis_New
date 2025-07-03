# Compiler and flags
CC = g++
CFLAGS = -std=c++17 -g -Wall

# Include paths
INCLUDES = -I $(PWD)/include \
           -I /u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/include/hipo4 \
           -I /u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.6/include \
           $(shell root-config --cflags)

# Linker flags
LDFLAGS = -L /u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/lib \
          -L /u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.6/lib \
          $(shell root-config --libs)

# Libraries to link against
LIBS = -lhipo4 -lccdb -pthread -lm -ldl

# Source and object files for skim
SKIM_SOURCES = skim.cpp
SKIM_OBJECTS = $(SKIM_SOURCES:.cpp=.o)

# Executable name
SKIM_EXEC = skim

# Default target
all: $(SKIM_EXEC)

# Link object files to create the skim executable
$(SKIM_EXEC): $(SKIM_OBJECTS)
	@echo "Linking objects for skim: $(SKIM_OBJECTS)"
	$(CC) $(SKIM_OBJECTS) -o $(SKIM_EXEC) $(LDFLAGS) $(LIBS)

# Compile source files to object files
%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Clean up
clean:
	rm -f $(SKIM_OBJECTS) $(SKIM_EXEC) *~ core

# Phony targets
.PHONY: all clean