# Pluto

[PLUTO][pluto] - An automatic parallelizer and locality optimizer for affine loop nests.

---

## Table of Contents

* [Installation](#installation)
    * [Prerequisites](#prerequisites)
    * [Building PLUTO](#building-pluto)
        * [Stable release](#stable-release)
        * [Development version from Git](#development-version-from-git)
    * [The polycc wrapper](#the-polycc-wrapper)
* [Command-line Options](#command-line-options)
* [Trying a new code](#trying-a-new-code)
* [Trying any included example](#trying-any-included-example)
* [More Information](#more-information)
* [Author](#author)
* [Contact](#contact)
* [License](#license)

---

## Installation

### Prerequisites

A Linux distribution. Pluto has been tested on x86 and x86-64 machines running
Fedora, Ubuntu, and RedHat Enterprise Server. Solaris should also be fine if you
have GNU utilities. 

In order to use the development version from Pluto's git repository, automatic build
system tools including *autoconf*, *automake*, and *libtool* are needed. *GMP* (GNU
multi precision arithmetic library) is needed by ISL (one of the included libraries).
If it's not already on your system, it can be installed easily with, for eg., 
`sudo yum -y install gmp gmp-devel` on a Fedora. 

It is also recommended `astyle` and `indent` be installed if a user wishes to browse 
through generated code.

Pluto includes all polyhedral libraries that it depends on.


### Building PLUTO

#### Stable release

```
tar zxvf pluto-0.11.4.tar.gz
cd pluto-0.11.4/
./configure
make
make test
```

The `configure` command can be provided `--with-isl-prefix=<isl install location>` to 
build with another isl, otherwise the bundled isl is used.

#### Development version from Git

```
git clone https://github.com/bondhugula/pluto.git
cd pluto/
git submodule init 
git submodule update
./autogen.sh
./configure [--enable-debug] [--with-isl-prefix=<isl install location>]
make
make test
```

The `configure` command can be provided `--with-isl-prefix=<location>` to compile
and link with an already installed isl. By default, the version of isl bundled
with Pluto will be used.

### The polycc wrapper

`polycc` is the wrapper script around src/pluto (core transformer) and all 
other components. `polycc` runs all of these in sequence on an input C 
program (with the section to  parallelize/optimize marked) and is what a 
user should use on input. Output generated is OpenMP parallel C code that 
can be readily compiled and run on shared-memory parallel machines like 
general-purpose multicores. `libpluto.{so,a}` is also built and can be found 
in `src/.libs/`. `make install` will install it.


## Command-line Options

To check the available command-line options please run:

```
./polycc -h
```

or see the documentation (`doc/DOC.txt`) for details.


## Trying a new code

- Use `#pragma scop` and `#pragma endscop` around the section of code 
  you want to parallelize/optimize.

- Then, just run 
    ```
    ./polycc <C source file> --parallel --tile
    ```

  The transformation is also printed out, and `test.par.c` will have the 
  parallelized code. If you want to see intermediate files, like the 
  .cloog file generated (`.opt.cloog`, `.tiled.cloog`, or `.par.cloog` 
  depending on command-line options provided), use `--debug` on command 
  line.

- Tile sizes can be specified in a file `tile.sizes`, otherwise default 
  sizes will be set. See `doc/DOC.txt` on how to specify the sizes.

To run a good number of experiments on a code, it is best to use the setup 
created for example codes in the `examples/` directory.  If you do not have 
ICC (Intel C compiler), uncomment line 7 and comment line 8 of 
`examples/common.mk` to use GCC.

- Just copy one of the sample directories in `examples/`, edit Makefile (SRC = )

- Do a make to build all executables
    - `orig` is the original code compiled with the native compiler
    - `tiled` is the tiled code
    - `par` is the OpenMP parallelized+locality optimized code
    - `lbpar` with diamond tiling  when possible
  One could do `make <target>` where target can be orig, orig_par, opt, tiled, par,
  lbpar, etc. (see `examples/common.mk` for full list)

- `make test` to test for correctness

- `make perf`, `make lbperf` to compare performance


## Trying any included example

Lets say we are trying the 2-d gauss seidel kernel. To access the application
please run:

```
$ cd examples/seidel
```

Inside the application's folder you can run several `make` commands to generate
the different codes. For instance:

- `make par` will generate `seidel.par.c` from `seidel.c` and also compile 
it to generate `par`.
- `make tiled` will generate `seidel.tiled.c` and also compile it to generate `tiled`
- `make orig` will compile `seidel.c` to generate `orig`.

Inside the application's folder you will find:

- `seidel.c`: This is the original code (the kernel in this code is extracted).  
'orig' is the corresponding executable when compiled with the native 
compiler (gcc or icc for eg.) with optimization flags, 'orig_par' with the 
native compiler's auto-parallelization enabled.
- `seidel.opt.c`: This is the transformed code without tiling (this is of not 
much use, except for seeing benefits of fusion in some cases). 'opt' is the 
corresponding executable.
- `seidel.tiled.c`: This is Pluto generated code optimized for locality with 
tiling and other transformations, but not not parallelized - this should be 
used for sequential execution. 'tiled' is the corresponding executable.
- `seidel.par.c`: This is Pluto parallelized code optimized for locality and 
parallelism  with tiling and other transformations. This code has OpenMP 
pragmas. 'par' is the corresponding executable.

To change any of the flags used for an example, edit the top section of 
`examples/common.mk` or the Makefile in the example directory

To manually specify tile sizes, create `tile.sizes`; see `examples/matmul/` 
for example or `doc/DOC.txt` for more information on setting tile sizes. 

The executables already have timers; you just have to run them and that will 
print execution time for the core part of the computation as well.

To run the Pluto parallelized version:
```
$ OMP_NUM_THREADS=4; ./par
```

To run native compiler optimized/auto-parallelized version:
```
$ OMP_NUM_THREADS=4; ./orig_par
```

To run the original unparallelized code:
```
$ ./orig
```

To run the locality optimized version generated by Pluto:
```
$ ./tiled
```

- `make clean` in the particular example's directory removes all executables 
    as well as generated codes

To launch a complete verification that compares output of tiled, par
with orig for all examples, in examples/, run `make test`.

```
[examples/ ]$ make test
```

## More Information

* See `doc/DOC.txt` for an overview of the system and details on all 
command-line options.
* For specifying custom tile sizes through `tile.sizes` file, see `doc/DOC.txt`.
* For specifying custom fusion structure through `.fst` file, see `doc/DOC.txt`.


## Author

Uday Bondhugula

uday@csa.iisc.ernet.in


## Contact

Please send all bug reports and comments to Uday Bondhugula 
<uday@csa.iisc.ernet.in> or post at pluto-development@googlegroups.com.

## License

Pluto is available under [GPL v3][gpl3], and libpluto is available under 
[LGPL v2.1][lgpl2.1].

[pluto]: http://pluto-compiler.sourceforge.net/
[gpl3]: https://www.gnu.org/licenses/gpl-3.0.en.html
[lgpl2.1]: https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
