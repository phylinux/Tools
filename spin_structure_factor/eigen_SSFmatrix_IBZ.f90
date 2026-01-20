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

	character(50)               :: filename, qlistfile
	character(100)              :: char1
	integer(4)                  :: Ntau
	complex(16), allocatable    :: SqSq(:,:,:,:,:,:)
	integer(4)                  :: qlist(1000, 3), nq
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
	integer(4)                  :: sb1p, sb2p
	real(8)                     :: clx(2)
	complex(8), allocatable     :: cmatrix(:,:)
	real(8)                     :: br
	complex(8), allocatable     :: diagval(:)
	real(8)                     :: aa(10,10)
	complex(8)                  :: zaa(10,10)
	real(8)                     :: eig(10)
	integer(4)                  :: i, j, k
	integer(4)                  :: iGN
	real(8)                     :: myzero(10)=0.d0

	integer(4)                  :: ixyz(3)
	real(8)                     :: Gvec(3)
	real(8)                     :: kvec(3), kvecp(3)
	real(8)                     :: kr

	narg = command_argument_count()
	if( narg /= 7 ) then
		write(*,*) "Please enter: filename, Ntau, Lx, Ly, Lz, sublat, qlistfile"
		!write(*,*) "Please enter: "
		!write(*,*) "case 1: filename, Ntau, Lx, Ly, Lz, sublat, nperiod"
		!write(*,*) "case 2: filename, Ntau, Lx, Ly, Lz, sublat, nperiod fm"
		!write(*,*) "case 3: filename, Ntau, Lx, Ly, Lz, sublat, nperiod, fm, h0, l0, m0"
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
	call getarg(7, qlistfile)

	Vol = Lx*Ly*Lz*sublat
	Lxyz = (/Lx, Ly, Lz/)
	allocate(SqSq(0:Ntau,sublat,sublat,0:Lx-1,0:Ly-1,0:Lz-1))
	SqSq = CMPLX(0.d0,0.d0,8)

	call read_ibz_info
	allocate(SqSq_IBZ(0:Ntau,sublat,sublat,ibz_size))
	SqSq_IBZ = CMPLX(0.d0,0.d0,8)

	allocate(sublatvec(sublat,3))
	latvec = 0.d0
	reclatvec = 0.d0
	sublatvec = 0.d0
	open(1, file="basis.bin", form="unformatted", status="old")
	read(1) latvec, sublatvec, reclatvec
	close(1)

	allocate(cmatrix(sublat,sublat))
	allocate(diagval(0:sublat))

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
	!open(21,file=trim(filename)//"-bin-trace.dat",action="write")
	!open(22,file=trim(filename)//"-bin-eigen.dat",action="write")
	do
		SqSq_IBZ = CMPLX(0.d0,0.d0)
		read(10,*,end=10) i
		!*** read one block ***!
		do i=1, ibz_size
			read(10,*) iq, ix, iy, iz
			dr = ibzw(i,1:3)-(/ix,iy,iz/)
			if( dot_product(dr,dr)/=0 ) then
				write(*,*) "Err: ibzw(i,1:3)-(/ix,iy,iz/) /= 0"
				stop
			end if
			do iGN=0, Ntau
				do sb2=1, sublat; do sb1=1,sublat      !!!!! sb2, sb1
					read(10,*) clx(1), clx(2)
					SqSq_IBZ(iGN,sb1,sb2,i) = CMPLX(clx(1),clx(2),8)
				end do; end do
			end do
		end do
		!*** one block done ***!

		do iq=1, nq
			rix = qlist(iq,1); ix = mod(rix,Lx)
			riy = qlist(iq,2); iy = mod(riy,Ly)
			riz = qlist(iq,3); iz = mod(riz,Lz)
			q = reclatvec(1,1:3)/Lx*rix + reclatvec(2,1:3)/Ly*riy + reclatvec(3,1:3)/Lz*riz
			write(20,'(4I4)') rix+riy*Lx+riz*Lx*Ly+1, rix, riy, riz
			!write(21,'(4I4)') rix+riy*Lx+riz*Lx*Ly+1, rix, riy, riz
			!write(22,'(4I4)') rix+riy*Lx+riz*Lx*Ly+1, rix, riy, riz
			!--- find the mapping point in IBZ ---!
			ixyz = (/ix,iy,iz/)  ! in BZ
			j = ix + iy*Lx + iz*Lx*Ly +1
			k = ibz_mapping(j)
			kvec  = 0.d0  ! k in BZ
			kvecp = 0.d0  ! k in IBZ
			Gvec  = 0.d0
			do i=1, 3
				kvec  = kvec  + reclatvec(i,1:3)/Lxyz(i)*ixyz(i)
				kvecp = kvecp + reclatvec(i,1:3)/Lxyz(i)*ibzw(k,i)
				Gvec  = Gvec  + reclatvec(i,1:3)        *ibz_G(j,i)
			end do
			!-------------------------------------!
			do iGN=0, Ntau
				diagval = CMPLX(0.d0, 0.d0)
				do sb1=1, sublat; do sb2=1, sublat
					sb1p = permutation(j,sb1)
					sb2p = permutation(j,sb2)

					br = dot_product(q,sublatvec(sb1,1:3)-sublatvec(sb2,1:3))
					br = -br
					clx(1) =  real(SqSq_IBZ(iGN,sb1p,sb2p,k))
					clx(2) = aimag(SqSq_IBZ(iGN,sb1p,sb2p,k))
					cmatrix(sb1,sb2) = CMPLX(clx(1),clx(2))*exp(CMPLX(0.d0,br))
					!--- extra phase due to rotation ---!
					kr = -dot_product(kvec,sublatvec(sb1,1:3)-sublatvec(sb2,1:3))
					cmatrix(sb1,sb2)  = cmatrix(sb1,sb2)/exp(CMPLX(0.d0,kr))
					kr =  dot_product(kvecp-Gvec,sublatvec(sb1p,1:3)-sublatvec(sb2p,1:3))
					cmatrix(sb1,sb2)  = cmatrix(sb1,sb2)/exp(CMPLX(0.d0,kr))
					!----------------------------------!
					diagval(0) = diagval(0) + cmatrix(sb1,sb2)
					!if( sb1==sb2 ) diagval(2) = diagval(2)+ cmatrix(sb1,sb2)
					if( sb1==sb2 ) diagval(sb1) = cmatrix(sb1,sb1)
				end do; end do

				!*** sum *******************************************************!
				write(20,'(10ES16.8)') real(diagval(0)), aimag(diagval(0)), myzero(1:10-2)
				!*** trace *****************************************************!
				!write(21,'(10ES16.8)') real(sum(diagval(1:sublat))), aimag(sum(diagval(1:sublat))), real(diagval(1:sublat)), myzero(1:10-sublat-2)
				!do sb2=1, sublat; do sb1=1, sublat
				!	write(21,'(2ES16.8)') cmatrix(sb1,sb2)
				!end do; end do
				!*** eigen *****************************************************!
				!=== real ==============================!
				aa = 0.d0
				!!--- do not use symmetry ---!
				!do sb1=1, sublat; do sb2=1, sublat
				!	aa(sb1,sb2) = real(cmatrix(sb1,sb2))
				!end do; end do
				!!--- use symmetry ---!
				!do sb1=1, sublat; do sb2=sb1, sublat
				!	aa(sb1,sb2) = real(cmatrix(sb1,sb2)+cmatrix(sb2,sb1))
				!	aa(sb2,sb1) = aa(sb1,sb2)
				!end do; end do
				!=== complex ===========================!
				zaa = CMPLX(0.d0,0.d0)
				!--- do not use symmetry ----!
				do sb1=1, sublat; do sb2=1, sublat
					zaa(sb1,sb2) = cmatrix(sb1,sb2)
				end do; end do
				!!--- use symmetry ----------!
				!do sb1=1, sublat; do sb2=sb1+1, sublat
				!	zaa(sb1,sb2) = (cmatrix(sb1,sb2)+conjg(cmatrix(sb2,sb1)))*0.5d0
				!	zaa(sb2,sb1) = conjg(zaa(sb1,sb2))
				!end do; end do
				!clx(1) = CMPLX(0.d0,0.d0)
				!do sb1=1, sublat
				!	clx(1) = clx(1) + zaa(sb1,sb1)
				!end do
				!clx(1) = clx(1)/sublat
				!do sb1=1, sublat
				!	zaa(sb1,sb1) = clx(1)
				!end do
				!write(*,*) "-------------------------------"
				!write(char1,'(I1,A25)') sublat,"(""("",F10.4,F10.4,"")"", 2X)"
				!do sb1=1, sublat
				!	write(*,'('//trim(char1)//')') cmatrix(sb1,1:sublat)
				!end do
				!write(*,*) "----"
				!do sb1=1, sublat
				!	write(*,'('//trim(char1)//')') zaa(sb1,1:sublat)
				!end do
				!write(*,*) "-------------------------------"
				!----------------!
				!call diasym(aa(1:sublat,1:sublat),eig(1:sublat),sublat)
				!call diazsym(zaa(1:sublat,1:sublat),eig(1:sublat),sublat)
				!call diaznsym(zaa(1:sublat,1:sublat),eig(1:sublat),sublat)
				!write(22,'(10ES16.8)') eig(1:sublat), sum(eig(1:sublat)), myzero(1:10-sublat-1)
				!write(*,'(8ES16.8)') eig(1:sublat), sum(eig(1:sublat)), myzero(1:8-sublat-1)
				!write(22,'(8ES16.8)') eig(1:sublat), sum(eig(1:sublat)), aimag(cmatrix(1,2)), aimag(cmatrix(1,3)), &
				!	& aimag(cmatrix(2,3)), myzero(1:8-sublat-4)
				!write(22,'(5ES16.8)') sum(eig(1:4)), eig(1:4)
			end do
		end do

		NBlck = NBlck+1
	end do
	10 close(10)
	close(20)
	!close(21)
	!close(22)
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
