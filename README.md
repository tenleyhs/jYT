Update 6/29/2018 by JLS

I am currently using this as a catch-all for my analysis scripts for simulations
run on hyades and fend. There are some duplications and there's some cleanup to
do, hopefully I'll get to this later.


Update 3/12/2012 by JFG

Brief instructions on installation:

* This package requires that YT is properly set up and enabled. Please read the
  YT installation guide for instructions on how to do this.

* It also requires that the "SciPy" package and the "progressbar" packages are
  installed. If the YT Python environment is enabled, these can be installed via
  the pip command. For SciPy, which depends on BLAS and LAPACK, you need to
  specify the location of the BLAS and LAPACK libraries before the pip command
  will work. Download the latest version of LAPACK (which includes BLAS), and
  your pip setup line should look something like this:

    CC=gcc FC=gfortran F77=gfortran BLAS_SRC=~/installers/lapack-3.4.0/BLAS/SRC/ LAPACK_SRC=~/installers/lapack-3.4.0/SRC/ pip install scipy

  It is important to select the compilers that were used to build YT on the
  system (which I think is gfortran by default).
