!--- infile is the output of SSF ---!

program structure_factor
	implicit none

	integer(4)                  :: narg                 ! the number of input argument

	character(50)               :: filename, qlistfile
	character(100)              :: char1
	integer(4)                  :: Ntau
	complex(16), allocatable    :: SqSq(:,:,:,:,:,:)
	integer(4)                  :: qlist(1000, 3), nq
	character(50)               :: charnum, charformat

	real(8)                     :: q(3)
	integer(4)                  :: Lx, Ly, Lz
	integer(4)                  :: subl
	integer(4)                  :: NBlck
	real(8)                     :: latvec(3,3)
	real(8)                     :: reclatvec(3,3)
	real(8), allocatable        :: sublatvec(:,:)

	integer(4)                  :: ix, iy, iz
	integer(4)                  :: rix, riy, riz
	integer(4)                  :: sb1, sb2
	real(8)                     :: clx(2)
	complex(8), allocatable     :: cmatrix(:,:)
	real(8)                     :: br
	!real(8)                     :: r(3)
	complex(8), allocatable     :: r(:)
	real(8)                     :: aa(10,10)
	complex(8)                  :: zaa(10,10)
	real(8)                     :: eig(10)
	integer(4)                  :: i, j
	integer(4)                  :: iGN
	real(8)                     :: myzero(10)=0.d0

	narg = command_argument_count()
	if( narg /= 7 ) then
		write(*,*) "Please enter: filename, Ntau, Lx, Ly, Lz, subl, qlistfile"
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
	call getarg(7, qlistfile)

	allocate(SqSq(0:Ntau,subl,subl,0:Lx-1,0:Ly-1,0:Lz-1))
	SqSq = CMPLX(0.d0,0.d0,8)

	allocate(cmatrix(subl,subl))
	allocate(r(0:subl))

	allocate(sublatvec(subl,3))
	latvec = 0.d0
	reclatvec = 0.d0
	sublatvec = 0.d0
	open(1, file="basis.bin", form="unformatted", status="old")
	read(1) latvec, sublatvec, reclatvec
	close(1)

	!--- Read qlist ------------------------------------------------------!
	open(1,file=trim(qlistfile),action="read")
	nq=0
	do
		read(1,*,end=1) qlist(nq+1,1:3)
		nq = nq+1
	end do
	1 close(1)
	!---------------------------------------------------------------------!

	NBlck = 0
	!--- START READ DATA -------------------------------------------------!
	open(10,file=trim(filename)//".txt",action="read")
	!open(10,file="../"//trim(filename)//".txt",action="read")
	!--- output bin ---!
	open(20,file=trim(filename)//"-bin-sum.dat",action="write")
	open(21,file=trim(filename)//"-bin-trace.dat",action="write")
	open(22,file=trim(filename)//"-bin-eigen.dat",action="write")
	do
		SqSq = CMPLX(0.d0,0.d0)
		read(10,*,end=10) i
		!*** read one block ***!
		do iz=0, Lz-1; do iy=0, Ly-1; do ix=0, Lx-1
			do iGN=0, Ntau
				do sb2=1, subl; do sb1=1,subl      !!!!! sb2, sb1
					read(10,*) clx(1), clx(2)
					SqSq(iGN,sb1,sb2,ix,iy,iz) = CMPLX(clx(1),clx(2),8)
				end do; end do
			end do
		end do; end do; end do
		!*** one block done ***!

		do j=1, nq
			rix = qlist(j,1); ix = mod(rix,Lx)
			riy = qlist(j,2); iy = mod(riy,Ly)
			riz = qlist(j,3); iz = mod(riz,Lz)
			q = reclatvec(1,1:3)/Lx*rix + reclatvec(2,1:3)/Ly*riy + reclatvec(3,1:3)/Lz*riz
			write(20,'(4I4)') rix+riy*Lx+riz*Lx*Ly+1, rix, riy, riz
			write(21,'(4I4)') rix+riy*Lx+riz*Lx*Ly+1, rix, riy, riz
			write(22,'(4I4)') rix+riy*Lx+riz*Lx*Ly+1, rix, riy, riz
			do iGN=0, Ntau
				r = CMPLX(0.d0, 0.d0)
				do sb1=1, subl; do sb2=1, subl
					br = dot_product(q,sublatvec(sb1,1:3)-sublatvec(sb2,1:3))
					br = -br
					clx(1) =  real(SqSq(iGN,sb1,sb2,ix,iy,iz))
					clx(2) = aimag(SqSq(iGN,sb1,sb2,ix,iy,iz))
					cmatrix(sb1,sb2) = CMPLX(clx(1),clx(2))*CMPLX(cos(br),sin(br))
					r(0) = r(0) + cmatrix(sb1,sb2)
					!if( sb1==sb2 ) r(2) = r(2)+ cmatrix(sb1,sb2)
					if( sb1==sb2 ) r(sb1) = cmatrix(sb1,sb1)
				end do; end do

				!*** sum *******************************************************!
				write(20,'(10ES16.8)') real(r(0)), aimag(r(0)), myzero(1:10-2)
				!*** trace *****************************************************!
				write(21,'(10ES16.8)') real(sum(r(1:subl))), aimag(sum(r(1:subl))), real(r(1:subl)), myzero(1:10-subl-2)
				!do sb2=1, subl; do sb1=1, subl
				!	write(21,'(2ES16.8)') cmatrix(sb1,sb2)
				!end do; end do
				!*** eigen *****************************************************!
				!=== real ==============================!
				aa = 0.d0
				!!--- do not use symmetry ---!
				!do sb1=1, subl; do sb2=1, subl
				!	aa(sb1,sb2) = real(cmatrix(sb1,sb2))
				!end do; end do
				!!--- use symmetry ---!
				!do sb1=1, subl; do sb2=sb1, subl
				!	aa(sb1,sb2) = real(cmatrix(sb1,sb2)+cmatrix(sb2,sb1))
				!	aa(sb2,sb1) = aa(sb1,sb2)
				!end do; end do
				!=== complex ===========================!
				zaa = CMPLX(0.d0,0.d0)
				!--- do not use symmetry ----!
				do sb1=1, subl; do sb2=1, subl
					zaa(sb1,sb2) = cmatrix(sb1,sb2)
				end do; end do
				!!--- use symmetry ----------!
				!do sb1=1, subl; do sb2=sb1+1, subl
				!	zaa(sb1,sb2) = (cmatrix(sb1,sb2)+conjg(cmatrix(sb2,sb1)))*0.5d0
				!	zaa(sb2,sb1) = conjg(zaa(sb1,sb2))
				!end do; end do
				!clx(1) = CMPLX(0.d0,0.d0)
				!do sb1=1, subl
				!	clx(1) = clx(1) + zaa(sb1,sb1)
				!end do
				!clx(1) = clx(1)/subl
				!do sb1=1, subl
				!	zaa(sb1,sb1) = clx(1)
				!end do
				!write(*,*) "-------------------------------"
				!write(char1,'(I1,A25)') subl,"(""("",F10.4,F10.4,"")"", 2X)"
				!do sb1=1, subl
				!	write(*,'('//trim(char1)//')') cmatrix(sb1,1:subl)
				!end do
				!write(*,*) "----"
				!do sb1=1, subl
				!	write(*,'('//trim(char1)//')') zaa(sb1,1:subl)
				!end do
				!write(*,*) "-------------------------------"
				!----------------!
				!call diasym(aa(1:subl,1:subl),eig(1:subl),subl)
				!call diazsym(zaa(1:subl,1:subl),eig(1:subl),subl)
				!call diaznsym(zaa(1:subl,1:subl),eig(1:subl),subl)
				!write(22,'(10ES16.8)') eig(1:subl), sum(eig(1:subl)), myzero(1:10-subl-1)
				!write(*,'(8ES16.8)') eig(1:subl), sum(eig(1:subl)), myzero(1:8-subl-1)
				!write(22,'(8ES16.8)') eig(1:subl), sum(eig(1:subl)), aimag(cmatrix(1,2)), aimag(cmatrix(1,3)), &
				!	& aimag(cmatrix(2,3)), myzero(1:8-subl-4)
				!write(22,'(5ES16.8)') sum(eig(1:4)), eig(1:4)
			end do
		end do

		NBlck = NBlck+1
	end do
	10 close(10)
	close(20)
	close(21)
	close(22)
	write(*,*) "NBlck = ", NBlck

END program structure_factor

!---------------------------!
subroutine diasym(aa,eig,n)
	!---------------------------!
	implicit none

	integer :: n,l,inf
	real(8) ::  aa(n,n),eig(n),work(n*(3+n/2))

	l=n*(3+n/2)
	call dsyev('V','U',n,aa,n,eig,work,l,inf)

end subroutine diasym
!---------------------!


!---------------------------!
!zheevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)
subroutine diazsym(aa,eig,n)
	!---------------------------!
	implicit none

	integer    :: n
	complex(8) :: aa(n,n)
	real(8)    :: eig(n)

	integer    ::  lwork
	integer    :: lrwork
	integer    :: liwork
	!complex(8)               ::   work(2*n+n**2)
	!real(8)                  ::  rwork(1+5*n+2*n**2)
	!integer                  ::  iwork(3+5*n)
	complex(8), allocatable  ::   work(:)
	real(8), allocatable     ::  rwork(:)
	integer, allocatable     ::  iwork(:)
	integer    :: info

    character(100) :: char1
	complex(8) :: aa_orig(n,n)
	complex(8) :: Av(n), lambda_v(n)
	real(8) :: error
	real(8), parameter :: tolerance = 1d-10
	integer(4)  :: i

	aa_orig = aa

    !write(*,*) "************* diazsym *************"
    !write(char1,'(I1,A25)') n,"(""("",F10.4,F10.4,"")"", 2X)"
    !do i=1, n
    !write(*,'('//trim(char1)//')') aa_orig(i,1:n)
    !end do

    != Workspace query 
	 lwork=-1
	lrwork=-1
	liwork=-1
	allocate(work(1))
	allocate(rwork(1))
	allocate(iwork(1))

	call zheevd('V','L',n,aa,n,eig,  work, lwork, &
	&                           rwork,lrwork, &
	&                           iwork,liwork, &
	& info)
	lwork = int(real(work(1)))
	lrwork = int(rwork(1))
	liwork = iwork(1)
	!write(*,*)  lwork, 2*n+n**2
	!write(*,*) lrwork, 1+5*n+2*n**2
	!write(*,*) liwork, 3+5*n
	!stop

	deallocate(work, rwork, iwork)
	allocate(work(lwork), rwork(lrwork), iwork(liwork))

	!! lwork=2*n+n**2
	!!lrwork=1+5*n+2*n**2
	!!liwork=3+5*n
	! lwork=198
	!lrwork=1+5*n+2*n**2
	!liwork=3+5*n
	!allocate(work(lwork), rwork(lrwork), iwork(liwork))

    aa = aa_orig
	call zheevd('V','L',n,aa,n,eig,  work, lwork, &
		&                           rwork,lrwork, &
		&                           iwork,liwork, &
		& info)

	!if (info /= 0) then
	!	print *, 'Error: ZHEEVD failed with INFO =', info
	!else
	!	print *, "Eigenvalues:"
	!	print '(6F12.6)', eig
	!	print *, "Eigenvectors (stored column-wise):"
	!	do i = 1, n
	!		print '(6("(",F8.4,",",F8.4,")"))', aa(1:n,i)
	!	end do
	!end if

	if( info/=0 ) then
		write(*,*) "Err: zheevd"
	else
		error = 0.d0
		do i = 1, n
			Av = matmul(aa_orig, aa(1:n,i))
			lambda_v = eig(i) * aa(1:n,i)
			error = error + maxval(abs(Av - lambda_v))
			if (error > tolerance) then
				write(*,'(6F12.8)') eig(1:6)
				write(*,*) "Err: zheevd, Av/=labmda*v", error
                stop
			end if
		end do
	end if


	deallocate(work, rwork, iwork)
end subroutine diazsym
!---------------------!


!---------------------------!
!zgeev(Jobvl, Jobvr, N, A, LDA, W, vl, ldvl, vr, ldvr, WORK, LWORK, RWORK, INFO)
subroutine diaznsym(aa,eig,n)
	!---------------------------!
	implicit none

	integer    :: n
	complex(8) :: aa(n,n)
	real(8)    :: eig(n)

	integer    :: lda
	integer    :: ldvl
	integer    :: ldvr
	integer    :: lwork
	complex(8), allocatable  ::   w(:)
	complex(8), allocatable  ::   work(:)
	real(8), allocatable     ::  rwork(:)
	complex(8), allocatable  ::   vl(:,:)
	complex(8), allocatable  ::   vr(:,:)
	integer    :: info

    character(100) :: char1
	complex(8) :: aa_orig(n,n)
	complex(8) :: Av(n), lambda_v(n)
	real(8) :: error
	real(8), parameter :: tolerance = 1d-10
	integer(4)  :: i

	aa_orig = aa
	lda = n
	ldvl = n
	ldvr = n

	allocate(w(n))
	allocate(vl(n,n))
	allocate(vr(n,n))
	allocate(rwork(2*n))

	!= Workspace query 
	lwork=-1
	allocate(work(1))

	call zgeev('N','V',n,aa,lda,w, vl,ldvl,vr,ldvr,  &
		& work, lwork, rwork, info)
	lwork = int(real(work(1))) +1
	!write(*,*)  lwork, 2*n+n**2
	!stop

	deallocate(work)
	allocate(work(lwork))

    aa = aa_orig
	call zgeev('N','V',n,aa,n,w, vl,ldvl,vr,ldvr,  &
		& work, lwork, rwork, info)
	eig(1:n) = real(w(1:n))

	if( info/=0 ) then
		write(*,*) "Err: zgeev"
	else
		error = 0.d0
		do i = 1, n
			Av = matmul(aa_orig, vr(1:n,i))
			lambda_v = w(i) * vr(1:n,i)
			error = error + maxval(abs(Av - lambda_v))
			if (error > tolerance) then
				write(*,'(6F12.8)') eig(1:6)
				write(*,*) "Err: zheevd, Av/=labmda*v", error
                stop
			end if
		end do
	end if


	deallocate(vl, vr, work, rwork)
end subroutine diaznsym
!---------------------!
