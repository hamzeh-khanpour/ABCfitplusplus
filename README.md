# ABCfit++
General constrained kinematic fit using ABC-parametrisation

## Installation

This is a shared library, in order to make the library you
can call
```
make
```

in the folder. This will create a library `libABCfit++.so`.

When in need of help execute `make help`.

## Usage

You can use the shared library by pointing in your gcc command lines
`-I` and `-L` to this directory or generally `-I` to these header
files and `-L` to the directory that the `libABCfit++.so` is
located. Then supply also `-lABCfit++` to the command line arguments.

An example:
```sh
g++ -o yourTestOutput -I../location/to/ABCfit++ -L../location/to/ABCfit++ -lABCfit++ YourSourceFile.cxx
```

## C++ versions
This library works with C++11 and above.

## Reference
This library is based on `ABCfit` developed by Oliver Buchmuller and Jørgen Beck Hansen. `ABCfit` is available [here](https://www.nbi.dk/~beck/abcfit/abcfit.html).
