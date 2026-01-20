!transform BZ data to IBZ data
!--- infile is the output of SSF ---!
!
! 0  !iblck
! iq ix iy iz
! real(SqSq(0,0,0))  aimag(SqSq(0,0,0))
! real(SqSq(1,0,0))  aimag(SqSq(1,0,0))
! real(SqSq(2,0,0))  aimag(SqSq(2,0,0))
! ...
! iq ix iy iz
! real(SqSq(0,0,0))  aimag(SqSq(0,0,0))
! real(SqSq(1,0,0))  aimag(SqSq(1,0,0))
! real(SqSq(2,0,0))  aimag(SqSq(2,0,0))
! ...
! 1  !iblck
! iq ix iy iz
! real(SqSq(0,0,0))  aimag(SqSq(0,0,0))
! real(SqSq(1,0,0))  aimag(SqSq(1,0,0))
! real(SqSq(2,0,0))  aimag(SqSq(2,0,0))
! ...
! iq ix iy iz
! real(SqSq(0,0,0))  aimag(SqSq(0,0,0))
! real(SqSq(1,0,0))  aimag(SqSq(1,0,0))
! real(SqSq(2,0,0))  aimag(SqSq(2,0,0))
! ...

program structure_factor
	implicit none

	integer(4)                  :: narg                 ! the number of input argument

	character(50)               :: filename
	character(100)              :: char1
	integer(4)                  :: Ntau
	complex(16), allocatable    :: SqSq(:,:,:,:,:,:)
	character(50)               :: charnum, charformat

	integer(4)                  :: ibz_size
	integer(4), allocatable     :: ibzw(:,:)
	integer(4), allocatable     :: ibz_mapping(:)
	integer(4), allocatable     :: ibz_G(:,:)
	real(8), allocatable        :: ibz_weight(:)
	integer(4), allocatable     :: permutation(:,:)
	complex(8), allocatable     :: SqSq_IBZ(:,:,:,:)

	real(8)                     :: q(3)
	integer(4)                  :: Vol
	integer(4)                  :: Lxyz(3)
	integer(4)                  :: Lx, Ly, Lz
	integer(4)                  :: sublat
	integer(4)                  :: NBlck
	real(8)                     :: latvec(3,3)
	real(8)                     :: reclatvec(3,3)
	real(8), allocatable        :: sublatvec(:,:)

	integer(4)                  :: iq
	integer(4)                  :: ix, iy, iz
	integer(4)                  :: rix, riy, riz
	integer(4)                  :: dr(3)
	integer(4)                  :: sb1, sb2
	real(8)                     :: clx(2)
	complex(8)                  :: c
	real(8)                     :: br
	integer(4)                  :: i, j, k
	integer(4)                  :: iGN

	integer(4)                  :: ixyz(3)
	real(8)                     :: Gvec(3)
	real(8)                     :: kvec(3), kvecp(3)
	real(8)                     :: kr
	integer(4)                  :: sb1p, sb2p

	narg = command_argument_count()
	if( narg /= 6 ) then
		write(*,*) "Please enter: filename, Ntau, Lx, Ly, Lz, sublat"
		stop
	end if

	call getarg(1, filename)
	call getarg(2, charnum)
	read(charnum,*) Ntau
	call getarg(3, charnum)
	read(charnum,*) Lx
	call getarg(4, charnum)
	read(charnum,*) Ly
	call getarg(5, charnum)
	read(charnum,*) Lz
	call getarg(6, charnum)
	read(charnum,*) sublat

	Vol = Lx*Ly*Lz*sublat
	Lxyz = (/Lx, Ly, Lz/)
	allocate(SqSq(0:Ntau,sublat,sublat,0:Lx-1,0:Ly-1,0:Lz-1))
	SqSq = CMPLX(0.d0,0.d0)

	call read_ibz_info
	allocate(SqSq_IBZ(0:Ntau,sublat,sublat,ibz_size))
	SqSq_IBZ = CMPLX(0.d0,0.d0)

	allocate(sublatvec(sublat,3))
	latvec = 0.d0
	reclatvec = 0.d0
	sublatvec = 0.d0
	open(1, file="basis.bin", form="unformatted", status="old")
	read(1) latvec, sublatvec, reclatvec
	close(1)


	NBlck = 0
	!--- START READ and TRANSFORM DATA ---------------------------------------!
	open(10,file=trim(filename)//".txt",action="read")
	open(20,file=trim(filename)//"-IBZ.txt",action="write")
	do
		SqSq = CMPLX(0.d0,0.d0)
		read(10,*,end=10) i
		!*** read one block ***!
		do iz=0, Lz-1; do iy=0, Ly-1; do ix=0, Lx-1
			do iGN=0, Ntau
				do sb2=1, sublat; do sb1=1,sublat      !!!!! sb2, sb1
					read(10,*) clx(1), clx(2)
					SqSq(iGN,sb1,sb2,ix,iy,iz) = CMPLX(clx(1),clx(2))
				end do; end do
			end do
		end do; end do; end do
		!*** one block done ***!

		SqSq_IBZ = CMPLX(0.d0,0.d0)
		!--- Irreducible Brillouin zone ---!
		do iz=0,Lz-1; do iy=0,Ly-1; do ix=0,Lx-1
			ixyz = (/ix,iy,iz/)
			j = ix + iy*Lx + iz*Lx*Ly +1
			k = ibz_mapping(j)
			kvec  = 0.d0
			kvecp = 0.d0
			Gvec  = 0.d0
			do i=1, 3
				kvec  = kvec  + reclatvec(i,1:3)/Lxyz(i)*ixyz(i)
				kvecp = kvecp + reclatvec(i,1:3)/Lxyz(i)*ibzw(k,i)
				Gvec  = Gvec  + reclatvec(i,1:3)        *ibz_G(j,i)
			end do
			do i=0, Ntau
				do sb2=1, sublat; do sb1=1, sublat
					sb1p = permutation(j,sb1)
					sb2p = permutation(j,sb2)
					if( sb1/=sb2 .and. sb1p==sb2p ) stop "Err: sb1/=sb2 .and. sb1p==sb2p"
					if( sb1==sb2 .and. sb1p/=sb2p ) stop "Err: sb1==sb2 .and. sb1p/=sb2p"
					c = SqSq(i,sb1,sb2,ix,iy,iz)*ibz_weight(j)
					kr = -dot_product(kvec,sublatvec(sb1,1:3)-sublatvec(sb2,1:3))
					c  = c*exp(CMPLX(0.d0,kr))
					kr =  dot_product(kvecp-Gvec,sublatvec(sb1p,1:3)-sublatvec(sb2p,1:3))
					c  = c*exp(CMPLX(0.d0,kr))
					SqSq_IBZ(i,sb1p,sb2p,k) = SqSq_IBZ(i,sb1p,sb2p,k) + c
				end do; end do
			end do
		end do; end do; end do

		NBlck = NBlck+1

		!- write to file -!
		write(20,'(I8)') NBlck
		do j = 1, ibz_size
			ix = ibzw(j,1)
			iy = ibzw(j,2)
			iz = ibzw(j,3)
			write(20,'(I8,3I4)') ix+iy*Lx+iz*Lx*Ly+1, ix, iy, iz
			do i=0, Ntau
				do sb2=1, sublat; do sb1=1, sublat
					c = SqSq_IBZ(i,sb1,sb2,j)
					write(20,'(2ES16.8)') real(c), aimag(c)
				end do; end do
			end do
		end do
	end do
	10 close(10)
	close(20)
	write(*,*) "NBlck = ", NBlck


CONTAINS

SUBROUTINE read_ibz_info
	implicit none
	integer(4)                  :: i, j
	integer(4)                  :: nr(3*3+1+sublat+1)

	open(1, file="ibzw.dat",action="read")
	read(1,*) ibz_size
	allocate(ibzw(ibz_size,4))
	do i=1, ibz_size
		read(1,*) ibzw(i,1:4)
	end do
	close(1)

	allocate(ibz_G(Vol/sublat,3))
	allocate(ibz_mapping(Vol/sublat))
	allocate(permutation(Vol/sublat,sublat))
	allocate(ibz_weight(Vol/sublat))
	open(1, file="ibz_mapping.dat",action="read")
	do i=1, Vol/sublat
		read(1,*) nr(1:3*3+1+sublat+1)
		j = nr(1) + nr(2)*Lx + nr(3)*Lx*Ly +1
		if( i/=j ) then
			write(*,*) "Err: ibz mapping, i/=j", i, j
			stop
		end if
		ibz_G(i,1:3) = nr(7:9)
		ibz_mapping(i) = nr(10)
		permutation(i,1:sublat) = nr(10+1:10+sublat)
		ibz_weight(i) = 1.d0/nr(10+sublat+1)
	end do
	close(1)
END SUBROUTINE read_ibz_info
END program structure_factor
