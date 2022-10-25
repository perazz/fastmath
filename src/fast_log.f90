!   ***********************************************************************************************
!   **                                                                                           **
!   **            |\   -  -.   ./                                                                **
!   **            | \./ \/ | ./ /     _________   _____________  ______  ________  __            **
!   **          __|         /  /     / ____/   | / ___/_  __/  |/  /   |/_  __/ / / /            **
!   **          \    .        /.-/  / /_  / /| | \__ \ / / / /|_/ / /| | / / / /_/ /             **
!   **           \   |\.|\/|    /  / __/ / ___ |___/ // / / /  / / ___ |/ / / __  /              **
!   **            \__\     /___/  /_/   /_/  |_/____//_/ /_/  /_/_/  |_/_/ /_/ /_/               **
!   **                                                                                           **
!   **                                         fastmath                                          **
!   **                                                                                           **
!   ***********************************************************************************************
!   **                                                                                           **
!   **    fast_log                                                                               **
!> @brief A module to compute FAST logarithm functions, based on Perini and Reitz, "Fast         **
!>        approximations of exponential and logarithm functions combined with efficient          **
!>        storage/retrieval for combustion kinetics calculations" Comb Flame 194(2018), 37-51.   **
!   **                                                                                           **
!   ***********************************************************************************************
!   **                                                                                           **
!> @author Copyright (C) Federico Perini <perini@wisc.edu>, 2016-2022
!> @since      tuesday,  10/01/2017
!   **                                                                                           **
!   ***********************************************************************************************
  module fast_log
     use iso_fortran_env, only: real64,int64,int16,int32,output_unit,real32
     implicit none
     private

     ! Module parameters
     integer, parameter      :: maxDegree    = 16
     real(real64), parameter :: half         = 0.5_real64
     real(real64), parameter :: one          = 1.0_real64
     real(real64), parameter :: two          = 2.0_real64
     real(real64), parameter :: log2         = log(two)
     real(real64), parameter :: rlog2        = one/log(two)
     real(real64), parameter :: sqrt2        = sqrt(two)

     integer(int64), parameter :: mantissa_left  = 2_int64**52
     integer(int64), parameter :: mantissa       = not(shiftl(2047_int64,52))
     integer(int64), parameter :: bias           = 1023_int64
     integer(int64), parameter :: ishift         = mantissa_left*bias

     public :: flog,flog_spline,flog_quintic
     public :: flog_name

     ! Module parameters
     integer, parameter, public :: FLOG_INTRINSIC  = 0
     integer, parameter, public :: FLOG_SPLINE3    = -1
     integer, parameter, public :: FLOG_SPLINE5    = -2

     ! A public type to compute fast exponentiatls
     type, public :: fast_logarithm
        integer :: type = FLOG_INTRINSIC
        contains

        procedure, non_overridable :: x     => flog_compute
        procedure, non_overridable :: set   => flog_set
        procedure, non_overridable :: print => flog_print

     end type fast_logarithm

     interface flog ! Generic has spline interpolation
        module procedure flog_spline
     end interface flog


     contains

     ! ********************************************************************************************
     ! TYPE-BOUND PROCEDURES
     ! ********************************************************************************************
     elemental real(real64) function flog_compute(this,x) result(log_x)
        class(fast_logarithm), intent(in) :: this
        real(real64), intent(in) :: x
        select case (this%type)
           case (FLOG_SPLINE3 ); log_x = flog_spline(x)
           case (FLOG_SPLINE5);  log_x = flog_quintic(x)
           case default;         log_x = log(x)
        end select
     end function flog_compute

     ! Set fast logarithm method
     elemental subroutine flog_set(this,degree)
        class(fast_logarithm), intent(inout) :: this
        integer, intent(in) :: degree
        select case (degree)
           case (-1);    this%type = FLOG_SPLINE3
           case (-2);    this%type = FLOG_SPLINE5
           case default; this%type = FLOG_INTRINSIC
        end select
     end subroutine flog_set

     ! Get name of the method
     pure function flog_name(flag)
        integer, intent(in) :: flag
        character(:), allocatable :: flog_name
        select case (flag)
           case (FLOG_SPLINE3); flog_name = 'spline(3)'
           case (FLOG_SPLINE5); flog_name = 'spline(5)'
           case default;        flog_name = 'intrinsic'
        end select
     end function flog_name

     ! Print fast logarithm details to external unit
     subroutine flog_print(this,iunit)
        class(fast_logarithm), intent(in) :: this
        integer, intent(in), optional :: iunit
        integer :: unitno
        if (present(iunit)) then
            unitno = iunit
        else
            unitno = output_unit
        endif

        write (unitno,'(A)') flog_name(this%type)

     end subroutine flog_print

     ! ********************************************************************************************
     ! DOUBLE PRECISION
     ! ********************************************************************************************
     elemental real(real64) function flog_spline(x) result(log_x)
        real(real64), intent(in) :: x

        real(real64) :: xi,xf
        integer(int64) :: i8
        real(real64), parameter :: s(3)= [rlog2,3.0_real64-2.5_real64*rlog2,1.5_real64*rlog2-2.0_real64]

        i8 = transfer(x,i8)
        xi = shiftr(i8,52)-bias

        ! Take mantissa part only
        xf = transfer(iand(i8,mantissa)+ishift,xf)-one

        ! Apply cubic polynomial
        xf = xf*(s(1)+xf*(s(2)+xf*s(3)))

        ! Compute log and Change of basis: log_2(x) -> log_e(x) = log2*log_2(x)
        log_x = (xf+xi)*log2

     end function flog_spline

     elemental real(real64) function flog_quintic(x) result(log_x)
        real(real64), intent(in) :: x

        real(real64) :: xi,xf
        integer(int64) :: i8
        real(real64), parameter :: s5(5)= [ 1.44269504088896e+0_real64,&
                                           -7.21347520444482e-1_real64,&
                                            4.42145354110618e-1_real64,&
                                           -2.12375830888126e-1_real64,&
                                            4.88829563330264e-2_real64]

        i8 = transfer(x,i8)
        xi = shiftr(i8,52)-bias

        ! Take mantissa part only
        xf = transfer(iand(i8,mantissa)+ishift,xf)-one

        ! Apply quintic polynomial
        xf = xf*(s5(1)+xf*(s5(2)+xf*(s5(3)+xf*(s5(4)+xf*s5(5)))))

        ! Compute log and Change of basis: log_2(x) -> log_e(x) = log2*log_2(x)
        log_x = (xf+xi)*log2

     end function flog_quintic

  end module fast_log

