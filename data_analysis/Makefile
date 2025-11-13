FCC=ifort
FCC=gfortran
FC=-O1 -cpp

F90=data_analysis.f90
EXEC=data_analysis

all:
	$(FCC) $(FC) data_analysis.f90      -o data_analysis
	$(FCC) $(FC) histogram_analysis.f90 -o histogram_analysis
	-rm -f *.mod
