! configuration in real space
! infile's format
!
! 0  !also could be other integer number
! Sx [Sy Sz]
! Sx [Sy Sz]
! ...
! 0
! Sx [Sy Sz]
! Sx [Sy Sz]
! ...

MODULE fftw_vrbls
	use, intrinsic :: iso_c_binding
	implicit none
#include <fftw3.f03>

	complex(C_DOUBLE_COMPLEX), allocatable :: fftw2_in(:,:), fftw2_out(:,:)
	complex(C_DOUBLE_COMPLEX), allocatable :: fftw3_in(:,:,:), fftw3_out(:,:,:)
	type(C_PTR)                            :: plan
END MODULE fftw_vrbls

PROGRAM r2q
	use fftw_vrbls
	implicit none

	integer(4)                  :: narg                 ! the number of input argument

	character(50)               :: filename
	integer(4)                  :: spintype
	real(8), allocatable        :: Sr(:,:,:,:,:)
	complex(8), allocatable     :: Sq(:,:,:,:,:)
	complex(8), allocatable     :: SqSq(:,:,:,:,:)
	integer(4)                  :: D
	integer(4)                  :: Lx, Ly, Lz
	integer(4)                  :: subl
	integer(4)                  :: NBlck
	real(8)                     :: reclatvec(3,3)
	real(8), allocatable        :: sublatvec(:,:)

	character(50)               :: charnum, charformat
	integer(4)                  :: ix, iy, iz
	integer(4)                  :: sb1, sb2
	real(8)                     :: c(3)
	integer(4)                  :: ndata
	integer(4)                  :: i, j

	narg = command_argument_count()
	if( narg/=7 ) then
		write(*,*) "Please enter: filename, spintype, D, Lx, Ly, Lz, subl"
		stop
	end if

	call getarg(1, filename)
	call getarg(2, charnum)
	read(charnum,*) spintype
	call getarg(3, charnum)
	read(charnum,*) D
	call getarg(4, charnum)
	read(charnum,*) Lx
	call getarg(5, charnum)
	read(charnum,*) Ly
	call getarg(6, charnum)
	read(charnum,*) Lz
	call getarg(7, charnum)
	read(charnum,*) subl

	if( D/=2 .and. D/=3 ) then
		write(*,*) "Err: D/= 2 or 3"
		stop
	end if

	call init_fftw

	allocate(Sr(subl, Lx,Ly,Lz,spintype))
	allocate(Sq(subl,Lx,Ly,Lz,spintype))
	allocate(SqSq(subl,subl,Lx,Ly,Lz))

	reclatvec = 0.d0
	open(10,file="reclatvec.dat",action="read")
	do i=1, 3
		read(10,*) reclatvec(i,1:3)
	end do
	close(10)

	allocate(sublatvec(subl,3))
	sublatvec = 0.d0
	open(10,file="sublatvec.dat",action="read")
	do i=1, subl
		read(10,*) sublatvec(i,1:3)
		sublatvec(i,1:3) = sublatvec(i,1:3) - (/1.d0,1.d0,1.d0/)*0.25d0
	end do
	close(10)

	NBlck = 0
	!SqSq = 0.d0
	SqSq = CMPLX(0.d0,0.d0)
	!--- START READ DATA -------------------------------------------------!
	open(10,file=trim(filename)//".txt",action="read")
	open(20,file="SqSq.txt",action="write")
	!- check data -!
	ndata = 0
	do
		read(10,*,end=10) i
		do iz=1, Lz; do iy=1, Ly; do ix=1, Lx
			do sb1=1, subl
				read(10,*) c(1:spintype)
				ndata = ndata+1
			end do
		end do; end do; end do
		NBlck = NBlck+1
	end do
	10 rewind(10)
	if( ndata/=(Lx*Ly*Lz*subl)*NBlck .or. &
		ix-1/=Lx .or. iy-1/=Ly .or. iz-1/=Lz .or. &
		sb1/=subl+1 ) then
		write(*,*) "Err: the data is incomplete."
		stop
	end if
	write(*,'(A10,I6)') "NBlck = ", NBlck
	!- read data -!
	NBlck = 0
	do
		read(10,*,end=20) i
		do iz=1, Lz; do iy=1, Ly; do ix=1, Lx
			do sb1=1,subl
				read(10,*) c(1:spintype)
				Sr(sb1, ix,iy,iz,1:spintype) = c(1:spintype)
			end do
		end do; end do; end do

		do i=1, spintype; do sb1=1,subl
			select case(D)
			case(2)
				do iy=1, Ly; do ix=1, Lx
					fftw2_in(ix,iy) = CMPLX(Sr(sb1,ix,iy,1,i),0.d0)
				end do; end do
				!-- FFTW --!
				call fftw_execute_dft(plan, fftw2_in, fftw2_out)
				Sq(sb1,:,:,1,i) = fftw2_out(:,:)
			case(3)
				do iz=1, Lz; do iy=1, Ly; do ix=1, Lx
					fftw3_in(ix,iy,iz) = CMPLX(Sr(sb1,ix,iy,iz,i),0.d0)
				end do; end do; end do
				!-- FFTW --!
				call fftw_execute_dft(plan, fftw3_in, fftw3_out)
				Sq(sb1,:,:,:,i) = fftw3_out(:,:,:)
			end select
		end do; end do

		do iz=1, Lz; do iy=1, Ly; do ix=1, Lx
			do sb1=1, subl;  do sb2=1, subl
				SqSq(sb1,sb2,ix,iy,iz) = dot_product(Sq(sb2,ix,iy,iz,1:spintype),Sq(sb1,ix,iy,iz,1:spintype))
			end do;  end do
		end do; end do; end do

		NBlck = NBlck+1
		write(20,'(I5)') NBlck
		do iz=1, Lz; do iy=1, Ly; do ix=1, Lx
			do sb2=1,subl; do sb1=1, subl
				write(20,'(2ES16.8)') Real(SqSq(sb1,sb2,ix,iy,iz)), Aimag(SqSq(sb1,sb2,ix,iy,iz))
			end do; end do
		end do; end do; end do
	end do

	20 close(10)
	close(20)
	!--- READ DATA END ---------------------------------------------------!

CONTAINS

SUBROUTINE init_fftw
	implicit none

	select case (D)
	case(2)
		call init_fftw2
	case(3)
		call init_fftw3
	case default
		write(*,*) "Err: D/=2,3"
		stop
	end select
END SUBROUTINE init_fftw

SUBROUTINE init_fftw2
	implicit none
	allocate(fftw2_in(Lx,Ly))
	allocate(fftw2_out(Lx,Ly))

	plan = fftw_plan_dft_2d(Ly,Lx, fftw2_in, fftw2_out, FFTW_FORWARD, FFTW_ESTIMATE)
END SUBROUTINE init_fftw2

SUBROUTINE init_fftw3
	implicit none
	allocate(fftw3_in(Lx,Ly,Lz))
	allocate(fftw3_out(Lx,Ly,Lz))

	plan = fftw_plan_dft_3d(Lz,Ly,Lx, fftw3_in, fftw3_out, FFTW_FORWARD, FFTW_ESTIMATE)
END SUBROUTINE init_fftw3
END PROGRAM r2q
