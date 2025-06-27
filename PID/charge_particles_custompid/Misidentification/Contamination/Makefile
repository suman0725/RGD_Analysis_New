# Compiler and flags
CC = g++
CFLAGS = -std=c++17 -g -Wall -Iinclude $(shell root-config --cflags)
LDFLAGS = $(shell root-config --libs)

# Source and output
SRC = main.C
EXEC = main

# Header files (in include/ directory)
HEADERS = include/BinUtils.h include/HistogramUtils.h include/FitUtils.h include/CanvasUtils.h include/LegendUtils.h include/DrawUtils.h include/ContaminationUtils.h

# Default target
all: $(EXEC)

# Build executable
$(EXEC): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SRC) -o $(EXEC) $(LDFLAGS)

# Clean up
clean:
	rm -f $(EXEC)

# Phony targets
.PHONY: all clean