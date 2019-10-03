gcc -lm -c lib_array.c lib_cpx.c lib_file.c lib_func.c -I../include
ar crv lib.a *.o
