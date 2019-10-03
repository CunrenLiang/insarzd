FFLAGS="-O2 -ffixed-line-length-132 -c"
gfortran $FFLAGS *.f
ar crv libinsarzdf.a *.o
ranlib libinsarzdf.a
