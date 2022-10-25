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
!   **    fast_rsqrt                                                                             **
!> @brief FAST reciprocal of a square root, 1/sqrt(x),  based on Perini and Reitz, "Fast         **
!>        approximations of exponential and logarithm functions combined with efficient          **
!>        storage/retrieval for combustion kinetics calculations" Comb Flame 194(2018), 37-51.   **
!   **                                                                                           **
!   ***********************************************************************************************
!   **                                                                                           **
!> @author Federico Perini <perini@wisc.edu>
!> @since      tuesday,  11/01/2017
!   **    Last modified: 11/01/2017                                                              **
!   **                                                                                           **
!   ***********************************************************************************************
  module fast_rsqrt
     use iso_fortran_env, only: real64,int64,int16,int32,output_unit,real32
     implicit none
     private

     ! Module parameters
     real(real64), parameter :: one          = 1.0_real64
     real(real64), parameter :: half         = 0.5_real64
     real(real64), parameter :: onep5        = 1.5_real64

     integer(int64), parameter :: magic      = 6910469410427058089_int64

     public :: frsqrt
     public :: frsqrt_name

     ! Module parameters
     integer, parameter, public :: FRSQRT_INTRINSIC = 0
     integer, parameter, public :: FRSQRT_QUAKE3    = 1

     ! A public type to compute fast exponentiatls
     type, public :: fast_invsqrt
        integer :: type = FRSQRT_INTRINSIC
        contains

        procedure, non_overridable :: x     => frsqrt_compute
        procedure, non_overridable :: set   => frsqrt_set
        procedure, non_overridable :: print => frsqrt_print

     end type fast_invsqrt

     interface frsqrt ! Generic has spline interpolation
        module procedure frsqrt_fast
     end interface frsqrt

     contains

     ! ********************************************************************************************
     ! TYPE-BOUND PROCEDURES
     ! ********************************************************************************************
     elemental real(real64) function frsqrt_compute(this,x) result(rsqrt_x)
        class(fast_invsqrt), intent(in) :: this
        real(real64), intent(in) :: x
        select case (this%type)
           case (FRSQRT_QUAKE3 ); rsqrt_x = frsqrt_fast(x)
           case default         ; rsqrt_x = one/sqrt(x)
        end select
     end function frsqrt_compute

     ! Set fast logarithm method
     elemental subroutine frsqrt_set(this,degree)
        class(fast_invsqrt), intent(inout) :: this
        integer, intent(in) :: degree
        select case (degree)
           case (1);     this%type = FRSQRT_QUAKE3
           case default; this%type = FRSQRT_INTRINSIC
        end select
     end subroutine frsqrt_set

     pure function frsqrt_name(flag)
        integer, intent(in) :: flag
        character(len=:), allocatable :: frsqrt_name

        select case (flag)
           case (FRSQRT_QUAKE3); frsqrt_name = 'Quake3'
           case default;         frsqrt_name = 'intrinsic'
        end select

     end function frsqrt_name

     ! Print fast logarithm details to external unit
     subroutine frsqrt_print(this,iunit)
        class(fast_invsqrt), intent(in) :: this
        integer, intent(in), optional :: iunit
        integer :: unitno
        if (present(iunit)) then
            unitno = iunit
        else
            unitno = output_unit
        endif

        write(unitno,*) frsqrt_name(this%type)

     end subroutine frsqrt_print

     ! ********************************************************************************************
     ! DOUBLE PRECISION
     ! ********************************************************************************************
     elemental real(real64) function frsqrt_fast(x) result(rsqrt_x)
        real(real64), intent(in) :: x

        real(real64) :: x2,y
        integer(int64) :: i

        x2 = half*x
        i  = transfer(x,i)
        i  = magic - shiftr(i,1)
        y  = transfer(i,rsqrt_x)

        ! Perform one Newton-Raphson step
        rsqrt_x  = y*(onep5-x2*y*y)

     end function frsqrt_fast

     ! ********************************************************************************************
     ! TEST
     ! ********************************************************************************************


  end module fast_rsqrt

