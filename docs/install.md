# Download and Installation

!!! abstract "Download the distribution"

    * [spgl1-2.1 (zip)](https://github.com/mpf/spgl1/archive/v2.1.zip)

## Installation

1. Unzip the distribution. This will create a directory called `spgl1-2.1/`.
2. Start Matlab and execute the following commands from the Matlab prompt.

```matlab
>> addpath spgl1-2.1  # Add Matlab to your path
>> cd spgl1-2.1       # Change directory
>> spgsetup           # Run SPGL1's setup routine
```

The spgsetup command compiles a fast C implementation of the projection
routines. Compiling Matlab MEX interfaces is sometimes tricky business, and if
your machine isn't setup for this, the `spgsetup` routine may fail. In that
case, SPGL1 will default to using the precompiled interfaces that have been
included. More information on how to change the default compiler is available
[here](https://www.mathworks.com/help/matlab/matlab_external/changing-default-compiler.html).



To verify that the SPGL1 installation is working, execute the following command from within Matlab:
```matlab
>> spgdemo
```

## Source code

The source code is maintained at https://github.com/mpf/spgl1.


