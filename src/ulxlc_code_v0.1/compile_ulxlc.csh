#!/bin/tcsh

ln -sf lmodel_ulxlc.dat lmodel.dat

set pkgname=ulxlc
echo "initpackage $pkgname lmodel.dat "`pwd`"\nquit\ny" | xspec
echo "lmod $pkgname .\nquit\ny" | xspec

rm -f *~ *.o *FunctionMap.* lpack_* *.mod Makefile
