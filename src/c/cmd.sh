#!/bin/bash


#files need library
SRCLIST1="extract_burst.c look.c look_msk.c mbf.c resamp.c rg_filter.c"
#files don't need library
SRCLIST2="correctphase.c flat.c interf.c mosaicframe.c mosaicsubswath.c psfilt1.c simamp.c"
#binary directory
BINDIR="../../bin"


if [ ! -d "$BINDIR" ]; then
    #echo "$BINDIR does not exist, create..."
    mkdir $BINDIR
fi


#FLAGS="-lm"
FLAGS="-D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -lm"
###################################################################
cd ./lib
./cmd.sh
cd ..

for srcf in ${SRCLIST1}
do
  #remove both directory and extension
  binf=`basename ${srcf} .c`
  gcc ${srcf} -I./include ./lib/lib.a -o ${BINDIR}/${binf} $FLAGS
done

cd ./lib
rm *.o *.a
cd ../

for srcf in ${SRCLIST2}
do
  #remove both directory and extension
  binf=`basename ${srcf} .c`
  gcc ${srcf} -o ${BINDIR}/${binf} $FLAGS
done

