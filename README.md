# fastmath

This is a Modern Fortran library for fast, approximate math functions: `exp`, `log`, `1/sqrt`. 

These functions provide fast, approximate evaluations of the exponential, logarithm and normalization functions in 64-bit precision. 

A detailed explanation of the underlying math with the IEEE-754 representation is reported in ref [1].

For each function, a module is provided; the interface has both simple functions and a programmable object-oriented wrapper. 

References
==========

This code is released under MIT license; please cite the following reference whenever using it.

[1] F. Perini, R.D. Reitz, "Fast approximations of exponential and logarithm functions combined with efficient storage/retrieval for combustion kinetics calculations" Combustion and Flame 194(2018), 37-51. 

Building, using
===============

An automated build is not available yet. 
- `src/fast_exp.f90`   contains the exponential routines
- `src/fast_log.f90`   contains the logarithm routines
- `src/fast_rsqrt.f90` contains the 1/sqrt(x) routines.
- `test/test.f90` is the test driver program 

Each module has no external dependencies, they can also be used individually.

A simple command line build script is: 

```
gfortran src/fast_exp.f90 src/fast_log.f90 src/fast_rsqrt.f90 test/test.f90  -o fastmath_test.exe
```

A simple makefile for the GNU compiler suite is provided in folder `project`; to run it: 

```
cd project/
make -f makefile.gcc
```

The testing executable compares performance and accuracy wrt the compiler's intrinsics functions. 
 
