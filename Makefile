CFLAGS=-fPIC -I include -Wall -pedantic
VER=c++11
CC=g++ -std=$(VER) $(CFLAGS)
LIBNAME=libABCfit++.so
SOURCES = $(wildcard src/*.cxx)
INCLUDES = $(wildcard include/*.h)
EXECUTABLES = $(patsubst examples/%.cxx,%,$(wildcard examples/*.cxx))
.PHONY: library examples help

library: build/$(LIBNAME)

help:
	@echo -----------------------------------------------------------------------
	@echo HELP
	@echo 
	@echo   library     : make shared library [default]
	@echo	examples    : compile and execute test cases from ABCclassesTest
	@echo   clean       : remove temporary files
	@echo -----------------------------------------------------------------------

build/%.o: src/%.cxx $(INCLUDES) | makethedir
	$(CC) -c -o $@ $< 

makethedir:
	@mkdir -p build

build/$(LIBNAME): $(patsubst src/%.cxx,build/%.o,$(SOURCES))
	$(CC) -shared -o $@ $^

%: examples/%.cxx build/$(LIBNAME) 
	$(CC) -o $@ $^ build/$(LIBNAME)

examples: $(EXECUTABLES)

clean:
	-rm build/$(LIBNAME)
	-rm build/*.o 
	-rm $(EXECUTABLES)

mrpropre: clean 
	-rm -r build