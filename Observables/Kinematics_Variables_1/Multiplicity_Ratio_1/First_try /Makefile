
CC=g++
CFLAGS=-std=c++11
INCLUDES=-I/home/suman/hipo/hipo4
LDFLAGS=-L/home/suman/hipo/lib -L/u/home/suman/hipo/lz4
LIBS=-lhipo4 -llz4 $(shell root-config --cflags --glibs)
%: %.cpp
	$(CC) $(CFLAGS) -o $@ $< $(INCLUDES) $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o *~ core
