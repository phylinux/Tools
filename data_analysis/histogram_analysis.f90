#include "Modules/stat_vrbls.f90"
#include "Modules/rng_vrbls.f90"

PROGRAM main
    use stat_vrbls
    use rng_vrbls
    implicit none

    integer(4)                  :: narg
    character(200)              :: filename
    character(100)              :: charmethod
    integer(4)                  :: method
    logical                     :: existyn
    integer(4)                  :: iblck
    integer(4)                  :: iobs

    real(8), allocatable        :: binObs(:,:)
    real(8), allocatable        :: initialdata(:,:,:)
    integer(4)                  :: len
    character(100)              :: charbins
    integer(4)                  :: bins

    real(8)                     :: maxv, minv, unitv
    real(8)                     :: c(2)
    integer(4)                  :: i, j

    narg = command_argument_count()

    select case (narg)
    case (1)
        call getarg(1,filename)
        method = 0
    case (3)
        call getarg(1,filename)
        call getarg(2,charmethod)
        read(charmethod,*) method
        call getarg(3,charbins)
        read(charbins,*) bins
    case default
        write(*,*) "Please enter the file name (and method, bins)."
        write(*,*) "OR if you want help, please use the parameter -h/--help"
        stop
    end select

    if( method/=0 .and. method/=1 ) then
        write(*,*) "method must be 0/1"
        stop
    end if

    inquire(file=trim(filename), exist=existyn)
    if( .not. existyn ) then
        write(*,*) "Err: No ", trim(filename)
        stop
    end if

    open(1,file=trim(filename), action="read")
    NBlck = 0
    read(1,*) len
    do
        read(1,*,end=10) c(1)
        do iobs=1, len
            read(1,*) c(1:2)
        end do
        NBlck = NBlck+1
    end do
    10 close(1)

    open(1,file="info.txt",action="write")
    write(1,'(A32,I6)') "The number of observables is", len
    write(1,'(A32,I6)') "The number of Block is", NBlck
    close(1)

    allocate(initialdata(len,NBlck,2))
    initialdata = 0.d0

    maxv = -1.d30
    minv =  1.d30
    open(1,file=trim(filename),action="read")
    read(1,*) len
    do iblck=1, NBlck
        read(1,*) c(1)
        do iobs=1, len
            read(1,*) c(1:2)
            initialdata(iobs,iblck,1:2) = c(1:2)
            if( maxv<c(1) ) maxv=c(1)
            if( minv>c(1) ) minv=c(1)
        end do
    end do
    close(1)

    NObs=500
    unitv = (maxv-minv)/NObs

    allocate(Obs(NObs,NBlck))
    Obs = 0.d0
    do iblck=1, NBlck
        do iobs=1, len
            c(1:2) = initialdata(iobs,iblck,1:2)
            j = int((c(1)-minv)/unitv)+1; if( J>NObs ) j=NObs
            Obs(j,iblck) = Obs(j,iblck) + c(2)
        end do
    end do
    !Obs = Obs*NObs

    if( method==1 ) then
        !write(*,'(A16,I8)') "Initial bins  = ", NBlck
        !write(*,'(A16,I8)') "Random  bins  = ", bins
        call date_and_time(date, time, zone, tval)
        Seed = tval(5)*3600+tval(6)*60+tval(7)
        call set_RNG

        allocate(binObs(NObs,bins))
        do iobs=1, NObs
            do iblck=1, bins
                i = int(rn()*NBlck)+1
                binObs(iobs,iblck) = Obs(iobs,i)
            end do
        end do
        NBlck = bins
        deallocate(Obs)
        allocate(Obs(NObs,NBlck))
        Obs = binObs
        deallocate(binObs)
    end if

    allocate(Ave(NObs))
    allocate(Dev(NObs))
    allocate(Cor(NObs))

    Ave = 0.d0
    Dev = 0.d0
    Cor = 0.d0

    call stat_analy(NBlck, method)

    do i = 1, NObs
    	write(*,'(4ES18.8)') minv+(i-1.d0+0.5d0)*unitv, Ave(i), Dev(i), Cor(i)
    enddo

    stop

CONTAINS
#include "Statistics/statistics.f90"
#include "WriteRead/wr_cnf_stat.f90"
#include "RNG/rng.f90"

END PROGRAM main
