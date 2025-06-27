# Compiler and flags
CC = g++
CFLAGS = -std=c++17 -g -Wall -Iinclude $(shell root-config --cflags)
LDFLAGS = $(shell root-config --libs)

# Source and output
SRC = main.cpp
EXEC = main

# Default target
all: $(EXEC)

# Build executable
$(EXEC): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(EXEC) $(LDFLAGS)

# Clean up
clean:
	rm -f $(EXEC)

# Phony targets
.PHONY: all clean