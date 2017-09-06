
LibForQ: A Fortran Library for Quantum Information Science

In this library we will continuously add and improve Fortran code to perform several numerical tasks one frequently 
needs when working in quantum information science. Among the functionalities already implemented, some examples are:
- Trace, partial trace, and partial transpose;
- Entanglement, discord, and coherence quantifiers;
- Pauli group (PG), Generalized Gell Mann Matrices (GGMM), and Bloch vector and correlation matrix with GGMM;
- Some other matrix functions are also available, some examples are:
  > Norms, inner products, and distance and distinguishability measures;
  > Purity, entropies, information measures, and mutual information;
  > Some popular, and not so popular, quantum states;
  > Array display, identity matrix, adjoint, Kronecker product (KP), KP of n elements of the PG,
    Gram-Schmidt orthogonalization, projector, outer product, etc.
I think it is better to take a look at the code and see all functions therein. The variables used and needed are 
explained there.

One can use this library simply by copying all the files to the main program folder and compiling them with:
$ gfortran -lblas -llapack *.f90
To run the executable file a.out generated just use:
$ ./a.out

Another, perhaps more convenient, manner of using the code is by creating a static LIBRARY. For that purpose you may 
follow the simple steps below:
1) Download the code
2) Go to the associated folder
3) Create a library with the commands:
$ gfortran -O3 -c *.f90
$ ar cr libforq.a *.o
To compile your main program using this library, copy libforq.a to your program's folder and use the command: 
$ gfortran -lblas -llapack libforq.a main.f90
Even better, you can also copy the library to your computer's libraries folder, e.g. with:
$ sudo cp libforq.a /usr/local/lib
and use it, anywhere in your computer, via
$ gfortran -lblas -llapack -lforq main.f90

REMARK #1: As seen above, this library depends on BLAS and Lapack. If you need random objects, 
see https://github.com/jonasmaziero/LibForro.
REMARK #2: Although on Mac all the commands above worked fine, on Ubuntu I had to use one of the following commands for compilation:
$ gfortran main.f90 /usr/local/lib/libforq.a /usr/lib/liblapack.a /usr/lib/libblas.a
$ gfortran main.f90 /usr/local/lib/libforq.a -llapack -lblas

Related references:
- J. Maziero, Random sampling of quantum states: A survey of methods, Braz. J. Phys. 45, 575 (2015), arXiv:1502.03644.
- J. Maziero, Generating pseudo-random discrete probability distributions,  Braz. J. Phys. 45, 377 (2015), arXiv:1502.02128.
- J. Maziero, Fortran code for generating random probability vectors, unitaries, and quantum states, Frontiers in ICT 3, 4 (2016)
  arXiv:1512.05173.
- J. Maziero, Computing partial traces and reduced density matrices, Int. J. Mod. Phys. C 28, 1750005 (2017). arXiv:1601.07458.
- J. Maziero, Computing coherence vectors and correlation matrices, with application to quantum discord quantification, 
  Adv. Math. Phys. 2016, 6892178 (2016), arXiv:1603.05284.
- J. Maziero, Computing partial transposes and related entanglement functions, Braz. J. Phys. 46, 605 (2016). arXiv:1609.00323.