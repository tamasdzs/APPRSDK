CC=g++
CFLAGS=-c -Wall -g -std=c++0x -I ../eigen -I. -I../cpp/src -I../src
LDFLAGS=
LIBS=../cpp/src/libalg.a
TESTS=approxTestWithHermite
SOURCES=$(TESTS:=.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=$(TESTS)

all: $(EXECUTABLE)

clean:
	rm -f $(EXECUTABLE) *.o
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
