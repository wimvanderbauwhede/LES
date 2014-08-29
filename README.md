# Large Eddy Simulator for the Study of Urban Boundary-layer Flows

This is an OpenCL version of the Large Eddy Simulator developed by Hiromasa Nakayama and  Haruyasu Nagai at the Japan Atomic Energy Agency and Prof. Tetsuya Takemi at the Disaster Prevention Research Institute of Kyoto University.

**NOTE:** This is _not_ the original Fortran-77 source code used in the included papers, but a Fortran-95 version where the complete time iteration loop is implemented in OpenCL. The papers are included to explain the physics behind the simulator.

## Prerequisites

To compile and run this code, you need to install the OclWrapper library from https://github.com/wimvanderbauwhede/OpenCLIntegration

You also need:

- a Fortran-95 compiler (tested with `gfortran` 4.8 and 4.9, `ifort` 12.0.0 and `pgf95` 12.5-0)
- a C++11-compliant C++ compiler for the OclWrapper (tested with `g++` 4.8 and 4.9)
- The `scons` build tool (I use v2.3) and therefore `python` v2.7
<!-- - Optionally, Perl if you want to generate a custom configuration -->
- Optionally, the NetCDF libraries from http://www.unidata.ucar.edu/software/netcdf/

## Compilation

Assuming you have installed the OclWrapper library correctly, compilation is simply:

    $ cd Sources
    $ scons

To compile the Fortran-95 code without the OclWrapper library,

    $ cd Sources
    $ scons -f SConstruct.F95_only

## Performance

I've tested this on a GeForce GTX 480 and a 12-core Intel Xeon, the OpenCL code is 7x faster on the GPU and 2x faster on the multicore CPU compared to the fastest compilation I could achieve, which was using the Intel Fortran compiler using auto-vectorisation and auto-parallelisation.

## Caveats

- At the moment the only way to configure the simulator is to edit the code.
- Also, the OpenCL code is entirely fixed, in terms of numbers of threads and compile time options. I will release the code generator for creating custom configurations soon.

## Copyright and License

This work is (c) Wim Vanderbauwhede, School of Computing Science, University of Glasgow
The source code may be used under the terms of either the

  * GNU Lesser General Public License (LGPL)
    http://www.opensource.org/licenses/lgpl-3.0.html

or the

  * BSD License (BSDL)
    http://www.opensource.org/licenses/bsd-license.php