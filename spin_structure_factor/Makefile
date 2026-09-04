
FCC=gfortran

F90=\
	SSF.f90

EXEC=SSF
FC=-O1 -cpp
FC=-O0 -g -cpp -DGFORT -fbounds-check -Wuninitialized

FFTW_INCLUDE_PATH=-I$(HOME)/opt/fftw-gcc/include
FFTW_LIBRARY_PATH=-L$(HOME)/opt/fftw-gcc/lib
FFTW_LIBS=-lfftw3

LAPACK_LIBRARY_PATH=-L$(HOME)/opt/lapack-gcc/lib
LAPACK_LIBS=-llapack -lblas

DLR_LIB_PATH=-L$(HOME)/opt/libdlr/lib
DLR_LIBS=-ldlr -Wl,-rpath,$(HOME)/opt/libdlr/lib

all: $(F90)
	$(FCC) $(FC) -o $(EXEC) $(F90)
	$(FCC) $(FC) SSF-sublattice.f90 -o SSF-sublattice
	#$(FCC) $(FC) SSF_IBZ.f90 -o SSF_IBZ
	$(FCC) $(FC) SSFmatrix_BZ2IBZ.f90 -o SSFmatrix_BZ2IBZ
	$(FCC) $(FC) r2q.f90 -o r2q $(FFTW_INCLUDE_PATH) $(FFTW_LIBRARY_PATH) $(FFTW_LIBS)
	$(FCC) $(FC) unform2form.f90 -o unform2form
	$(FCC) $(FC) eigen_SSFmatrix.f90 $(LAPACK_LIBS) -o eigen_SSFmatrix
	$(FCC) $(FC) eigen_SSFmatrix_IBZ.f90 $(LAPACK_LIBS) -o eigen_SSFmatrix_IBZ
	$(FCC) $(FC) dsyev.f XXYY.f90 -o XXYY
	$(FCC) $(FC) eigen_SSFmatrix_IBZ_DLR.f90 $(DLR_LIB_PATH) $(DLR_LIBS) -o eigen_SSFmatrix_IBZ_DLR

clean:
	-rm -rf $(EXEC) *.mod
