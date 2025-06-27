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

# ROOT-only libraries for analyzeContamination
ROOT_LIBS = $(shell root-config --libs) -pthread -lm -ldl

# Source and object files for skim
SKIM_SOURCES = skim.cpp src/ParticleData.cpp
SKIM_OBJECTS = $(SKIM_SOURCES:.cpp=.o)

# Source and object files for plot_chi2pid
PLOT_SOURCES = plot_chi2pid.cpp
PLOT_OBJECTS = $(PLOT_SOURCES:.cpp=.o)

# Source and object files for chi2pid_310_345
CHI2PID_SOURCES = chi2pid_310_345.cpp
CHI2PID_OBJECTS = $(CHI2PID_SOURCES:.cpp=.o)

# Source and object files for test
TEST_SOURCES = test.cpp src/ParticleData.cpp
TEST_OBJECTS = $(TEST_SOURCES:.cpp=.o)

# Source and object files for analyzeContamination
CONTAMINATION_SOURCES = analyzeContamination.cpp src/BinUtils.cc src/HistogramUtils.cc \
                        src/CanvasManager.cc src/FitFunctions.cc src/ContaminationUtils.cc
CONTAMINATION_OBJECTS = $(patsubst %.cpp, %.o, $(filter %.cpp, $(CONTAMINATION_SOURCES))) \
                       $(patsubst src/%.cc, src/%.o, $(filter %.cc, $(CONTAMINATION_SOURCES)))

# Executable names
SKIM_EXEC = skim
PLOT_EXEC = plot_chi2pid
CHI2PID_EXEC = chi2pid_310_345
TEST_EXEC = test
CONTAMINATION_EXEC = analyzeContamination

# Default target
all: $(SKIM_EXEC) $(PLOT_EXEC) $(CHI2PID_EXEC) $(TEST_EXEC) $(CONTAMINATION_EXEC)

# Link object files to create the skim executable
$(SKIM_EXEC): $(SKIM_OBJECTS)
	@echo "Linking objects for skim: $(SKIM_OBJECTS)"
	$(CC) $(SKIM_OBJECTS) -o $(SKIM_EXEC) $(LDFLAGS) $(LIBS)

# Link object files to create the plot_chi2pid executable
$(PLOT_EXEC): $(PLOT_OBJECTS)
	@echo "Linking object for plot_chi2pid: $(PLOT_OBJECTS)"
	$(CC) $(PLOT_OBJECTS) -o $(PLOT_EXEC) $(LDFLAGS) $(LIBS)

# Link object files to create the chi2pid_310_345 executable
$(CHI2PID_EXEC): $(CHI2PID_OBJECTS)
	@echo "Linking object for chi2pid_310_345: $(CHI2PID_OBJECTS)"
	$(CC) $(CHI2PID_OBJECTS) -o $(CHI2PID_EXEC) $(LDFLAGS) $(LIBS)

# Link object files to create the test executable
$(TEST_EXEC): $(TEST_OBJECTS)
	@echo "Linking object for test: $(TEST_OBJECTS)"
	$(CC) $(TEST_OBJECTS) -o $(TEST_EXEC) $(LDFLAGS) $(LIBS)

# Link object files to create the analyzeContamination executable
$(CONTAMINATION_EXEC): $(CONTAMINATION_OBJECTS)
	@echo "Linking objects for analyzeContamination: $(CONTAMINATION_OBJECTS)"
	$(CC) $(CONTAMINATION_OBJECTS) -o $(CONTAMINATION_EXEC) $(ROOT_LIBS)

# Compile source files to object files
%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

src/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

src/%.o: src/%.cc
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Dependencies
analyzeContamination.o: analyzeContamination.cpp include/CanvasManager.h include/BinUtils.h include/HistogramUtils.h include/FitFunctions.h include/ContaminationUtils.h
src/CanvasManager.o: src/CanvasManager.cc include/CanvasManager.h include/BinUtils.h
src/BinUtils.o: src/BinUtils.cc include/BinUtils.h
src/HistogramUtils.o: src/HistogramUtils.cc include/HistogramUtils.h include/BinUtils.h
src/FitFunctions.o: src/FitFunctions.cc include/FitFunctions.h
src/ContaminationUtils.o: src/ContaminationUtils.cc include/ContaminationUtils.h

# Clean up
clean:
	rm -f $(SKIM_OBJECTS) $(PLOT_OBJECTS) $(CHI2PID_OBJECTS) $(TEST_OBJECTS) $(CONTAMINATION_OBJECTS) \
	      $(SKIM_EXEC) $(PLOT_EXEC) $(CHI2PID_EXEC) $(TEST_EXEC) $(CONTAMINATION_EXEC) *~ core

# Phony targets
.PHONY: all clean