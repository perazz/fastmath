program test_fastmath
    use fast_exp
    use fast_log
    use fast_rsqrt
    use iso_fortran_env
    implicit none

    integer :: nfailed=0,npassed=0

    call add_test(test_fast_exp())
    call add_test(test_fast_log())
    call add_test(test_fast_rsqrt())

    print 1, npassed+nfailed,npassed,nfailed
    if (nfailed>0) then
        stop -1
    else
        stop 0
    endif

    1 format('[fastmath] ',i0,' test completed: ',i0,' passed, ',i0,' failed.')

    contains

    subroutine add_test(success)
        logical, intent(in) :: success
        if (success) then
            npassed = npassed+1
        else
            nfailed = nfailed+1
        end if
    end subroutine add_test

    ! Test fast 1/sqrt(x)
    logical function test_fast_rsqrt() result(success)

        integer, parameter :: nsize = 1000000
        integer, parameter :: ntest = 500
        real(real64), parameter :: lxmin =  -200_real64
        real(real64), parameter :: lxmax =  +200_real64
        real(real64), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i,method
        real(real64) :: time,c_start,c_end
        real(real64), dimension(FRSQRT_INTRINSIC:FRSQRT_QUAKE3) :: timep,avgerr

        allocate(x(nsize),intrin(nsize),packge(nsize),z(ntest))

        call random_number(x)

        ! Do a logspace
        x = 10.0_real64**(lxmin*(1.0_real64-x) + lxmax*x)

        ! Get intrinsic
        time = 0.0_real64
        do i=1,ntest
            call cpu_time(c_start)
            intrin = 1.0_real64/sqrt(x)
            call cpu_time(c_end)
            z(i) = sum(intrin)
            time = time+c_end-c_start
        end do

        print *, '*** fast 1/sqrt(x) test ***'

        do method=FRSQRT_INTRINSIC,FRSQRT_QUAKE3

            ! Set method
            timep(method) = 0.0_real64
            do i=1,ntest
                call cpu_time(c_start)
                select case (method)
                   case (FRSQRT_QUAKE3);    packge = frsqrt(x)
                   case (FRSQRT_INTRINSIC); packge = 1.0_real64/sqrt(x)
                end select
                call cpu_time(c_end)
                z(i) = sum(packge)
                timep(method) = timep(method)+c_end-c_start
            end do

            ! Get average error
            avgerr(method) = sum(abs(packge-intrin)/abs(intrin),intrin/=0.0_real64)/count(intrin/=0.0_real64)

            ! Print stats
            print 1, frsqrt_name(method), &
                     1e9*timep(method)/(nsize*ntest),&
                     time/timep(method),&
                     avgerr(method)

        end do

        success = maxval(avgerr,1)<=1.0e-2_real64

        1 format(a24,': average time = ',f9.4,' ns/eval, speed-up=',f5.2,'X, relative error=',1pe15.5e3)

    end function test_fast_rsqrt


    ! Test fast logarithm
    logical function test_fast_log() result(success)

        integer, parameter :: nsize = 100000
        integer, parameter :: ntest = 200
        real(real64), parameter :: lxmin =  -200_real64
        real(real64), parameter :: lxmax =  +200_real64
        real(real64), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i,method
        real(real64) :: time,c_start,c_end
        real(real64), dimension(FLOG_SPLINE5:FLOG_INTRINSIC) :: timep,avgerr

        allocate(x(nsize),intrin(nsize),packge(nsize),z(ntest))

        call random_number(x)

        ! Do a logspace
        x = 10.0_real64**(lxmin*(1.0_real64-x) + lxmax*x)

        ! Get intrinsic
        time = 0.0_real64
        do i=1,ntest
            call cpu_time(c_start)
            intrin = log(x)
            call cpu_time(c_end)
            z(i) = sum(intrin)
            time = time+c_end-c_start
        end do

        print *, '*** fast log(x) test ***'

        do method=FLOG_SPLINE5,FLOG_INTRINSIC

            ! Set method
            timep(method) = 0.0_real64
            do i=1,ntest
                call cpu_time(c_start)
                select case (method)
                   case (FLOG_SPLINE3); packge = flog_spline(x)
                   case (FLOG_SPLINE5); packge = flog_quintic(x)
                   case default;        packge = log(x)
                end select
                call cpu_time(c_end)
                z(i) = sum(packge)
                timep(method) = timep(method)+c_end-c_start
            end do

            ! Get average error
            avgerr(method) = sum(abs(packge-intrin)/abs(intrin),intrin/=0.0_real64)/count(intrin/=0.0_real64)

            ! Print stats
            print 1, flog_name(method), &
                     1e9*timep(method)/(nsize*ntest),&
                     time/timep(method),&
                     avgerr(method)

        end do

        success = maxval(avgerr,1)<=1.0e-4_real64

        1 format(a24,': average time = ',f9.4,' ns/eval, speed-up=',f5.2,'X, relative error=',1pe15.5e3)

    end function test_fast_log

    ! Test fast exponential
    logical function test_fast_exp() result(success)

        integer, parameter :: nsize = 200000
        integer, parameter :: ntest = 300
        real(real64), parameter :: xmin = -300.0_real64
        real(real64), parameter :: xmax =  300.0_real64
        real(real64), allocatable :: x(:),intrin(:),packge(:),z(:)
        integer :: i,method
        real(real64) :: time,c_start,c_end
        real(real64), dimension(FEXP_SPLINE5:FEXP_POLY5) :: timep,avgerr

        allocate(x(nsize),intrin(nsize),packge(nsize),z(ntest))

        call random_number(x)
        x = xmin*(1.0_real64-x) + xmax*x

        ! Get intrinsic
        time = 0.0_real64
        do i=1,ntest
            call cpu_time(c_start)
            intrin = exp(x)
            call cpu_time(c_end)
            z(i) = sum(intrin)
            time = time+c_end-c_start
        end do

        print *, '*** fast exp(x) test ***'

        do method=FEXP_SPLINE5,FEXP_POLY5

            ! Set method
            timep(method) = 0.0_real64
            do i=1,ntest
                call cpu_time(c_start)
                select case (method)
                   case (FEXP_POLY1  ); packge = fexp1(x)
                   case (FEXP_POLY2  ); packge = fexp2(x)
                   case (FEXP_POLY3  ); packge = fexp3(x)
                   case (FEXP_POLY4  ); packge = fexp4(x)
                   case (FEXP_POLY5  ); packge = fexp5(x)
                   case (FEXP_SPLINE3); packge = fexps3(x)
                   case (FEXP_SPLINE5); packge = fexps5(x)
                   case default       ; packge = exp(x)
                end select
                call cpu_time(c_end)
                z(i) = sum(packge)
                timep(method) = timep(method)+c_end-c_start
            end do

            ! Get average error
            avgerr(method) = sum(abs(packge-intrin)/abs(intrin),intrin/=0.0_real64)/count(intrin/=0.0_real64)

            ! Print stats
            print 1, fexp_name(method), &
                     1e9*timep(method)/(nsize*ntest),&
                     time/timep(method),&
                     avgerr(method)

        end do

        success = maxval(avgerr,1)<=5.0e-2_real64

        1 format(a24,': average time = ',f9.4,' ns/eval, speed-up=',f5.2,'X, relative error=',1pe15.5e3)

    end function test_fast_exp

end program test_fastmath
