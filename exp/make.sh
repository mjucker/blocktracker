#!/bin/bash

rm *.x *.o *.mod

export FLAGS='-g -traceback -r8 -Duse_netCDF -I/sayre/include/ -L/sayre/lib -lnetcdff -lnetcdf'

ifort $FLAGS -c ../src/blob_tracker.f90
ifort $FLAGS -c ../src/conrec.f90
ifort $FLAGS -c ../example/main.f90

ifort $FLAGS conrec.o blob_tracker.o main.o -o ex.x


