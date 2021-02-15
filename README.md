# ljmd

Contributors:

 - Avinash
anand-avinash

-Giulio Dondi
giuliodondi

- Matteo Poggi
pogmat


Compilation of LJMD code:
the basic compilation steps are:

- mkdir -p build
- cd build
- cmake (specify_flags_here) ..
- cmake --build .

The flags that may be specified are:

- -DCMAKE_BUILD_TYPE=type : will label the build type of the project
- -DENABLE_TIMING=ON : will enable more detialed sampling of the run-times of various portions of the code: i/o read, integration, force calculation, kinetic energy calculation.  
In the MPI version the times sampled by all the processes will be averaged at the end.
- -DENABLE_TESTING=ON : compiles unit tests with the googletest library

- -DENABLE_MPI=ON : will compile the MPI version
- -DENABLE_OMP=ON : will compile the OMP version

In the current version only one between OMP and MPI should be chosen.
