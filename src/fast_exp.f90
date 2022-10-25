!   ***********************************************************************************************
!   **                                                                                           **
!   **                  |\   -  -.   ./                                                          **
!   **                  | \./ \/ | ./ /     __________  ___________ __________                   **
!   **                __|         /  /     / ____/ __ \/ ____/ ___// ____/ __ \                  **
!   **                \    .        /.-/  / /_  / /_/ / __/  \__ \/ /   / / / /                  **
!   **                 \   |\.|\/|    /  / __/ / _, _/ /___ ___/ / /___/ /_/ /                   **
!   **                  \__\     /___/  /_/   /_/ |_/_____//____/\____/\____/                    **
!   **                                                                                           **
!   **                                     FRESCO-SpeedCHEM                                      **
!   **            A code for internal combustion engine flows with chemical kinetics             **
!   **                                                                                           **
!   ***********************************************************************************************
!   **                                                                                           **
!   **    fast_exp                                                                               **
!> @brief A module to compute FAST exponential functions, based on Perini and Reitz, "Fast       **
!>        approximations of exponential and logarithm functions combined with efficient          **
!>        storage/retrieval for combustion kinetics calculations" Comb Flame 194(2018), 37-51.   **
!   **                                                                                           **
!   ***********************************************************************************************
!   **                                                                                           **
!> @author Federico Perini <perini@wisc.edu>
!> @since      monday ,  27/06/2016
!   **    Last modified: 27/06/2016                                                              **
!   **                                                                                           **
!   ***********************************************************************************************
  module fast_exp
     use iso_fortran_env, only: real64,int64,output_unit
     implicit none
     private

     ! Module parameters
     integer :: ideg
     integer, parameter      :: maxDegree    = 16
     real(real64), parameter :: half         = 0.5_real64
     real(real64), parameter :: one          = 1.0_real64
     real(real64), parameter :: two          = 2.0_real64
     real(real64), parameter :: log2         = log(two)
     real(real64), parameter :: rlog2        = one/log(two)
     real(real64), parameter :: sqrt2        = sqrt(two)

     ! Series centered in x=1/2
     real(real64), parameter :: x0 = half
     real(real64), parameter :: f0 = 1.5_real64-sqrt2
     real(real64), parameter :: k(maxDegree) = [one-sqrt2*log2,&
                                                (-sqrt2*log2**ideg/gamma(ideg+one),ideg=2,maxDegree)]
     integer(int64), parameter :: mantissa = 2_int64**52
     integer(int64), parameter :: bias     = 1023_int64
     integer(int64), parameter :: ishift   = mantissa*bias

     public :: fexp,fexp1,fexp2,fexp3,fexp4,fexp5,fexps3,fexps5
     public :: fexp_name

     ! Module parameters
     integer, parameter, public :: FEXP_INTRINSIC  =  0
     integer, parameter, public :: FEXP_SPLINE3    = -1
     integer, parameter, public :: FEXP_SPLINE5    = -2
     integer, parameter, public :: FEXP_POLY1      =  1
     integer, parameter, public :: FEXP_POLY2      =  2
     integer, parameter, public :: FEXP_POLY3      =  3
     integer, parameter, public :: FEXP_POLY4      =  4
     integer, parameter, public :: FEXP_POLY5      =  5

     ! A programmable object-oriented wrapper
     type, public :: fast_exponential
        integer   :: type = FEXP_SPLINE3
        contains

        procedure, nopass :: x_spline => fexp_spline
        procedure, nopass :: x_quintic => fexp_quintic

        procedure, non_overridable :: x      => fexp_compute
        procedure, non_overridable :: set    => fexp_set
        procedure, non_overridable :: print  => fexp_print

     end type fast_exponential

     interface fexp ! Generic has 5th-order interpolation
        module procedure fexp5_double
     end interface fexp
     interface fexp1
        module procedure fexp1_double
     end interface fexp1
     interface fexp2
        module procedure fexp2_double
     end interface fexp2
     interface fexp3
        module procedure fexp3_double
     end interface fexp3
     interface fexp4
        module procedure fexp4_double
     end interface fexp4
     interface fexp5
        module procedure fexp5_double
     end interface fexp5
     interface fexpn
        module procedure fexpn_double
     end interface fexpn
     interface fexps3
        module procedure fexp_spline
     end interface fexps3
     interface fexps5
        module procedure fexp_quintic
     end interface fexps5

     contains

     ! ********************************************************************************************
     ! TYPE-BOUND PROCEDURES
     ! ********************************************************************************************
     elemental real(real64) function fexp_compute(this,x) result(exp_x)
        class(fast_exponential), intent(in) :: this
        real(real64), intent(in) :: x
        select case (this%type)
           case (FEXP_POLY1  ); exp_x = fexp1(x)
           case (FEXP_POLY2  ); exp_x = fexp2(x)
           case (FEXP_POLY3  ); exp_x = fexp3(x)
           case (FEXP_POLY4  ); exp_x = fexp4(x)
           case (FEXP_POLY5  ); exp_x = fexp5(x)
           case (FEXP_SPLINE3); exp_x = fexp_spline(x)
           case (FEXP_SPLINE5); exp_x = fexp_quintic(x)
           case default            ; exp_x = exp(x)
        end select
     end function fexp_compute

     ! Set fast exponential method
     elemental subroutine fexp_set(this,degree)
        class(fast_exponential), intent(inout) :: this
        integer, intent(in) :: degree
        select case (degree)
           case (1:5);   this%type = degree
           case (-1);    this%type = FEXP_SPLINE3
           case (-2);    this%type = FEXP_SPLINE5
           case default; this%type = FEXP_INTRINSIC
        end select
     end subroutine fexp_set

     ! Get method name
     pure function fexp_name(flag)
        integer, intent(in) :: flag
        character(len=:), allocatable :: fexp_name

        select case (flag)
           case (FEXP_POLY1);    fexp_name = 'linear   accuracy'
           case (FEXP_POLY2);    fexp_name = 'degree 2 accuracy'
           case (FEXP_POLY3);    fexp_name = 'degree 3 accuracy'
           case (FEXP_POLY4);    fexp_name = 'degree 4 accuracy'
           case (FEXP_POLY5);    fexp_name = 'degree 5 accuracy'
           case (FEXP_SPLINE3);  fexp_name = 'cubic spline accuracy'
           case (FEXP_SPLINE5);  fexp_name = 'quintic spline accuracy'
           case default;         fexp_name = 'intrinsic'
        end select

     end function fexp_name

     ! Print fast exponential details to external unit
     subroutine fexp_print(this,iunit)
        class(fast_exponential), intent(in) :: this
        integer, intent(in), optional :: iunit
        integer :: unitno
        if (present(iunit)) then
            unitno = iunit
        else
            unitno = output_unit
        endif

        write(unitno,'(A)')fexp_name(this%type)

     end subroutine fexp_print

     ! ********************************************************************************************
     ! DOUBLE PRECISION
     ! ********************************************************************************************

     elemental real(real64) function fexp5_double(x) result(exp_x)
        real(real64), intent(in) :: x
        exp_x = fexpn_double(x,5)
     end function fexp5_double

     elemental real(real64) function fexp4_double(x) result(exp_x)
        real(real64), intent(in) :: x
        exp_x = fexpn_double(x,4)
     end function fexp4_double

     elemental real(real64) function fexp3_double(x) result(exp_x)
        real(real64), intent(in) :: x
        exp_x = fexpn_double(x,3)
     end function fexp3_double

     elemental real(real64) function fexp2_double(x) result(exp_x)
        real(real64), intent(in) :: x
        exp_x = fexpn_double(x,2)
     end function fexp2_double

     elemental real(real64) function fexp1_double(x) result(exp_x)
        real(real64), intent(in) :: x
        exp_x = fexpn_double(x,1)
     end function fexp1_double

     elemental real(real64) function fexpn_double(x,n) result(exp_x)
        real(real64), intent(in) :: x
        integer     , intent(in) :: n

        real(real64) :: xr,xf,xi
        integer(int64) :: i8
        integer :: i

        xr = x*rlog2 ! Change of basis: 2^x -> e^x
        xf = xr-floor(xr)

        xi = xf-x0
        do i=1,min(n,maxDegree)
          xr = xr - xi*k(i)
          xi = xi*(xf-x0)
        end do
        xr = xr - f0

        i8 = int(mantissa*xr + ishift,int64)
        exp_x = transfer(i8,exp_x)

     end function fexpn_double

     ! Use a 3rd-degree spline, differentiable at the extremes
     elemental real(real64) function fexp_spline(x) result(exp_x)
        real(real64), intent(in) :: x

        real(real64) :: xr,xf
        integer(int64) :: i8

        ! Cubic spline expansion coefficients
        real(real64), parameter :: s(3)= [-7.6167541742324804e-2_real64,&
                                          -0.23141283591588344_real64  ,&
                                           0.30758037765820823_real64   ]

        xr = x*rlog2         ! Change of base: e^x -> 2^y
        xf = xr-floor(xr)   ! Compute fractional part
        xr = xr-((s(1)*xf+s(2))*xf+s(3))*xf ! Subtract Delta
        i8 = int(mantissa*xr + ishift,int64) !
        exp_x = transfer(i8,exp_x)

     end function fexp_spline

     ! Use a 5th-degree spline, twice differentiable at the extremes
     elemental real(real64) function fexp_quintic(x) result(exp_x)
        real(real64), intent(in) :: x

        real(real64) :: xr,xf
        integer(int64) :: i8

        real(real64), parameter :: s5(5)= [-1.90188191959304e-3_real64,&
                                           -9.01146535969578e-3_real64,&
                                           -5.57129652016652e-2_real64,&
                                           -2.40226506959101e-1_real64,&
                                            3.06852819440055e-1_real64]

        xr = x*rlog2 ! Change of base: e^x -> 2^y
        xf = xr-floor(xr)
        xr = xr-((((s5(1)*xf+s5(2))*xf+s5(3))*xf+s5(4))*xf+s5(5))*xf
        i8 = int(mantissa*xr + ishift,int64)
        exp_x = transfer(i8,exp_x)

     end function fexp_quintic

  end module fast_exp
