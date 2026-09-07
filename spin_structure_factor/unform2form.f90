
program main
	implicit none
	integer(4)                  :: narg                 ! the number of input argument

	character(50)               :: filename
	complex(8), allocatable     :: SxSx(:,:,:,:,:,:)
	integer(4)                  :: Ntau
	integer(4)                  :: Lx, Ly, Lz
	integer(4)                  :: subl

	character(50)               :: charnum
	integer(4)                  :: iblck
	complex(8)                  :: c
	integer(4)                  :: sb1, sb2
	integer(4)                  :: ix, iy, iz
	integer(4)                  :: i

	narg = command_argument_count()
	if( narg<6 ) then
		write(*,*) "Please enter: filename, Ntau, Lx, Ly, Lz, subl"
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
	read(charnum,*) subl


	allocate(SxSx(0:Ntau,subl,subl,0:Lx-1,0:Ly-1,0:Lz-1))

	open(1,file=trim(filename)//".txt",form="UNFORMATTED",access="sequential",status="old",action="read")
	open(2,file=trim(filename)//"-formatted.txt",action="write")
	iblck=0
	do
		read(1,end=10) SxSx
		iblck = iblck+1

		write(2,'(I8)') iblck
		do iz=0,Lz-1; do iy=0,Ly-1; do ix=0,Lx-1
			!write(1,'(I8,3I4)') ix+iy*Lx+iz*Lx*Ly+1, ix, iy, iz
			do i=0, Ntau
				do sb2=1, subl; do sb1=1, subl
					c = SxSx(i,sb1,sb2,ix,iy,iz)
					write(2,'(2ES16.8)') real(c), aimag(c)
				end do; end do
			end do
		end do; end do; end do
	end do
	10 close(1)
	close(2)
	write(*,*) "NBlck = ", iblck

	deallocate(SxSx)

end program main
