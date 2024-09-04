
FCC=gfortran

F90=\
	SSF.f90

EXEC=SSF
FC=-O1 -cpp
#FC=-O0 -g -cpp -DGFORT -fbounds-check -Wuninitialized

MY_INCLUDE_PATH=-I$(HOME)/opt/fftw-gcc/include
MY_LIBRARY_PATH=-L$(HOME)/opt/fftw-gcc/lib
LIBS=-lfftw3

all: $(F90)
	$(FCC) $(FC) -o $(EXEC) $(F90)
	$(FCC) $(FC) r2q.f90 -o r2q $(MY_INCLUDE_PATH) $(MY_LIBRARY_PATH) $(LIBS)

clean:
	-rm -rf $(EXEC) *.mod
