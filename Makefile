CFLAGS=-fPIC
VER=c++11
CC=g++ -std=$(VER) $(CFLAGS)
LIBNAME=libABCfit++.so

.PHONY: library test help

library: $(LIBNAME)

help:
	@echo -----------------------------------------------------------------------
	@echo HELP
	@echo 
	@echo   library : make shared library [default]
	@echo	test    : compile and execute test cases from ABCclassesTest
	@echo   clean   : remove temporary files
	@echo -----------------------------------------------------------------------

%.o: %.cxx %.h
	$(CC) -c -o $@ $<

$(LIBNAME): $(patsubst %.cxx,%.o,$(wildcard *.cxx))
	$(CC) -shared -o $@ $^

test: $(patsubst %.cxx,%.o,$(wildcard *.cxx))
	mkdir -p test/bin/
	$(CC) -o test/bin/test $^
	test/bin/./test

clean:
	rm -f $(LIBNAME)
	rm -f *.o
	rm -rf test
