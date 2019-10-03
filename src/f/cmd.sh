#!/bin/bash


#files need library
SRCLIST1="fitoff.f get_peg.f rect.f rect_with_looks.f"
#files don't need library
SRCLIST2=""
#binary directory
BINDIR="../../bin"


if [ ! -d "$BINDIR" ]; then
    #echo "$BINDIR does not exist, create..."
    mkdir $BINDIR
fi


#
FFLAGS="-O2 -ffixed-line-length-132"
###################################################################
cd ./libf
./cmd.sh
cd ../

for srcf in ${SRCLIST1}
do
  #remove both directory and extension
  binf=`basename ${srcf} .f`
  gfortran $FFLAGS -c -o ${binf}.o ${srcf}
  gfortran $FFLAGS -o ${BINDIR}/${binf} ${binf}.o ./libf/libinsarzdf.a
  rm ${binf}.o
done

cd ./libf
rm *.o *.a
cd ../

for srcf in ${SRCLIST2}
do
  #remove both directory and extension
  binf=`basename ${srcf} .f`
  gfortran $FFLAGS -o ${BINDIR}/${binf} ${srcf}
done

