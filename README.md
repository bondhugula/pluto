# Pluto

## OVERVIEW

Please see http://pluto-compiler.sourceforge.net.

This package includes both the tool pluto and libpluto. The `pluto` tool is a
source-to-source transformer meant to be run via the polycc script, `libpluto`
provides a thread-safe library interface.

[![Pluto build and
test](https://github.com/bondhugula/pluto/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/bondhugula/pluto/actions/workflows/build_and_test.yml)

[![Check format with clang-format](https://github.com/bondhugula/pluto/actions/workflows/clang_format.yml/badge.svg)](https://github.com/bondhugula/pluto/actions/workflows/clang_format.yml)

## LICENSE

Pluto and libpluto are available under the MIT LICENSE. Please see the file
`LICENSE` in the top-level directory for more details.

## INSTALLING PLUTO

### PREREQUISITES

A Linux distribution. Pluto has been tested on x86 and x86-64 machines running
Fedora, Ubuntu, and CentOS.

- In order to use the development version from Pluto's git repository, automatic
build system tools, including `autoconf`, `automake`, and `libtool` are needed.

- LLVM/Clang 14.x (14.x recommended, 11.x, 12.x tested to work as well), along
  with its development/header files, is needed for the pet submodule. These
  packages are available in standard distribution repositories or could be
  installed by building LLVM and Clang from sources. See `pet/README` for
  additional detail.  On most modern distributions, these can be installed from
  the repositories.

  Example:
  ```shell
  # On an Ubuntu.
  sudo apt install -y llvm-14-dev libclang-14-dev
  # On a Fedora.
  sudo dnf -y install llvm14-devel clang14-devel
  ```

- LLVM `FileCheck` is used for Pluto's test suite. (On a Fedora, this is part of
  the 'llvm' package.)

- GMP (GNU multi precision arithmetic library) is needed by ISL (one of the
  included libraries).  If it's not already on your system, it can be installed
  easily with, for eg., `sudo yum -y install gmp gmp-devel` on a Fedora (`sudo
  apt-get install libgmp3-dev` or something similar on an Ubuntu).

Pluto includes all polyhedral libraries that it depends on. See `pet/README` for
pet's pre-requisites.

### BUILDING PLUTO

**Stable release:**

```shell
$ tar zxvf pluto-0.11.4.tar.gz
$ cd pluto-0.11.4/
$ ./configure [--with-clang-prefix=<clang install location>]
$ make
$ make test
```

configure can be provided `--with-isl-prefix=<isl install location>` to build
with another isl, otherwise the bundled isl is used.

**Development version from Git:**

```shell
git clone git@github.com:bondhugula/pluto.git
cd pluto/
git submodule init
git submodule update
./autogen.sh
./configure [--enable-debug] [--with-clang-prefix=<clang headers/libs location>]
# Example: on an Ubuntu: --with-clang-prefix=/usr/lib/llvm-14, on a Fedora,
# typically, it's /usr/lib64/llvm14.
make
make check-pluto
```

* Use `--with-clang-prefix=<location>` to point to the specific clang to
build with.

* Use `--with-isl-prefix=<isl install location>` to compile and link with an
already installed isl. By default, the version of isl bundled with Pluto will be
used.

`polycc` is the wrapper script around src/pluto (core transformer) and all other
components. `polycc` runs all of these in sequence on an input C program (with
the section to parallelize/optimize marked) and is what a user should use on
input. Output generated is OpenMP parallel C code that can be readily compiled
and run on shared-memory parallel machines like general-purpose multicores.
`libpluto.{so,a}` is also built and can be found in `src/.libs/`. `make install`
will install it.

## TRYING A NEW EXAMPLE

- Use `#pragma scop` and `#pragma endscop` around the section of code
  you want to parallelize/optimize.

- Then, just run `./polycc <C source file>`.

  The transformation is also printed out, and `test.par.c` will have the
  parallelized code. If you want to see intermediate files, like the
  `.cloog` file generated (`.opt.cloog`, `.tiled.cloog`, or `.par.cloog`
  depending on command-line options provided), use `--debug` on command
  line.

- Tile sizes can be specified in a file `tile.sizes`, otherwise default
  sizes will be set. See `doc/DOC.txt` on how to specify the sizes.

To run a good number of experiments on a code, it is best to use the setup
created for example codes in the `examples/` directory.  If you do not have
`ICC` (Intel C compiler), uncomment line 9 and comment line
8 of `examples/common.mk` to use GCC.

- Just copy one of the sample directories in `examples/`, edit `Makefile` (`SRC
= `).

- do a `make` (this will build all executables; `orig` is the original code
compiled with the native compiler, `tiled` is the tiled code, `par` is the
OpenMP parallelized + locality optimized code. One could do `make <target>`
where target can be orig, orig_par, opt, tiled, par, pipepar, etc (see
`examples/common.mk` for full list).

- `make check-pluto` to test for correctness, `make perf` to compare
performance.

## COMMAND-LINE OPTIONS

Run

```shell
./polycc -h
```

or see documentation (`doc/DOC.txt`) for details.


## TRYING ANY INCLUDED EXAMPLE CODE

Let's say we are trying the 2-d gauss seidel kernel. In `examples/seidel`, do
`make par`; this will generate `seidel.par.c` from `seidel.c` and also compile
it to generate `par`.  Likewise, `make tiled` for `tiled` and `make orig` for
`orig`.

```shell
cd examples/seidel
```

`seidel.c`: This is the original code (the kernel in this code is extracted).
`orig` is the corresponding executable when compiled with the native compiler
(`gcc` or `icc` for eg.) with optimization flags, `orig_par` with the native
compiler's auto-parallelization enabled.

`seidel.opt.c`: This is the transformed code without tiling (this is of not much
use, except for seeing benefits of fusion in some cases). `opt` is the
corresponding executable.

`seidel.tiled.c`: This is Pluto generated code optimized for locality with
tiling and other transformations, but not not parallelized - this should be used
for sequential execution. `tiled` is the corresponding executable.

`seidel.par.c`: This is Pluto parallelized code optimized for locality and
parallelism  with tiling and other transformations. This code has OpenMP
pragmas. `par` is the corresponding executable.

- To change any of the flags used for an example, edit the top section of
`examples/common.mk` or the `Makefile` in the example directory

- To manually specify tile sizes, create `tile.sizes`; see `examples/matmul/`
for example or `doc/DOC.txt` for more information on setting tile sizes.

The executables already have timers; you just have to run them and that will
print execution time for the core part of the computation as well.

To run the Pluto parallelized version:

```shell
OMP_NUM_THREADS=4 ./par
```

To run native compiler optimized/auto-parallelized version:

```shell
OMP_NUM_THREADS=4 ./orig_par
```

To run the original unparallelized code:

```shell
./orig
```

To run the locality optimized version generated by Pluto:

```shell
./tiled
```

`make clean` in the particular example's directory removes all executables as
well as generated codes.

To launch a complete verification that compares output of tiled, par with orig
for all examples, in `examples/`, run `make check-pluto`.

```shell
[examples/ ]$ make check-pluto
```

## MORE INFORMATION

* See `doc/DOC.txt` for an overview of the system and details on all
command-line options.

## BUGS AND ISSUES

Please report bugs and issues at https://github.com/bondhugula/pluto/issues.

For questions and general discussion, please email
pluto-development@googlegroups.com after joining the group:
https://groups.google.com/g/pluto-development.
