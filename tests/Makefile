CC=g++
CFLAGS=-c -Wall -g -std=c++0x -I ../eigen -I. -I../cpp/src -I../src -I ../plotLib/matplotlib-cpp -I/usr/include/python2.7
LDFLAGS=
LIBS=../cpp/src/libalg.a -lpython2.7
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
