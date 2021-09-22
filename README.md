# ABCfit++
General constrained kinematic fit using ABC-parametrisation

## Installation

This is a shared library, in order to make the library you can call
```
make
```

in the folder. This will create a library `libABCfit++.so` in the (automatic created) build folder.

When in need of help execute `make help`.


## Usage

The project comes with two demo/test examples to demonstrate the usage. The command
```
make examples
```
will create the executables `demo` and `ABCtestclasses`. The former can be consulted for inspiration.

To employ in your own code, you can use the shared library by pointing in your g++ command lines
`-I` and `-L` to this directory or generally `-I` to these header
files and `-L` to the directory that the `libABCfit++.so` is
located. Then supply also `-lABCfit++` to the command line arguments.

An example:
```sh
g++ -std=c++11 -fPIC -I location/to/ABCfitplusplus/include -L location/to/ABCfitplusplus/build -lABCfit++ -o YourExeFile  YourSourceFile.cxx
```

Depending on your OS you might need to update the environment variable `LD_LIBRAY_PATH`to include the location of the libABCfit++.so.

## C++ versions
This library works with C++11 and above.

## Reference
This library is based on `ABCfit` developed by Oliver Buchmuller and Jørgen Beck Hansen. `ABCfit` is available [here](https://www.nbi.dk/~beck/abcfit/abcfit.html).
