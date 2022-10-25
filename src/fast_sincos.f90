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
!   **    fast_exp                                                                               **
!> @brief Fast sin/cos calculations: work in progress
!>        This package is based on Perini and Reitz, "Fast approximations of exponential and     **
!>        logarithm functions combined with efficient storage/retrieval for combustion kinetics  **
!>        calculations" Comb Flame 194(2018), 37-51.                                             **
!   **                                                                                           **
!   ***********************************************************************************************
!   **                                                                                           **
!> @author Copyright (C) Federico Perini <perini@wisc.edu>, 2016-2022
!> @since      monday ,  27/06/2016
!   **                                                                                           **
!   ***********************************************************************************************
module fast_sincos
    use iso_fortran_env
    implicit none

    public :: sincos32

    contains

     elemental real(real32) function cos_32s(x) result(cosx)
        real(real32), intent(in) :: x
        real(real32), parameter :: c1= 0.99940307
        real(real32), parameter :: c2=-0.49558072
        real(real32), parameter :: c3= 0.03679168
        real(real32) :: x2      ! The input argument squared
        x2=x*x
        cosx = c1 + x2*(c2 + c3 * x2)
     end function cos_32s

     elemental real(real32) function fast_cos(angle) result(cosx)
       real(real32), intent(in) :: angle

       real(real32), parameter :: invtwopi=0.1591549
       real(real32), parameter :: pi      =3.141593
       real(real32), parameter :: twopi   =6.283185
       real(real32), parameter :: halfpi  =1.570796
       real(real32), parameter :: threehalfpi=4.7123889
       real(real32) :: angle2pi

       ! clamp to the range 0..2pi
       angle2pi=angle-floor(angle*invtwopi)*twopi
       angle2pi=merge(angle2pi,-angle2pi,angle2pi>0.0)

       if(angle2pi<halfpi) then
          cosx = cos_32s(angle2pi)
       elseif (angle2pi<pi) then
          cosx = cos_32s(pi-angle2pi)
       elseif (angle2pi<threehalfpi) then
          cosx = cos_32s(angle2pi-pi)
       else
          cosx = cos_32s(twopi-angle2pi)
       end if
     end function fast_cos

     ! Fast sin+cos calculation in single precision
     elemental subroutine sincos32(x,sinx,cosx)
        real(real32), intent(in) :: x
        real(real32), intent(out) :: sinx,cosx

        real(real32), parameter :: invtwopi   =0.1591549
        real(real32), parameter :: pi         =3.141593
        real(real32), parameter :: twopi      =6.283185
        real(real32), parameter :: halfpi     =1.570796
        real(real32), parameter :: threehalfpi=4.7123889
        real(real32) :: angle2pi,sinmultiplier

        ! clamp to the range 0..2pi
        angle2pi     =x-floor(x*invtwopi)*twopi
        sinmultiplier=merge(1.0,-1.0,angle2pi>0.0 .and. angle2pi<pi)
        angle2pi     =merge(angle2pi,-angle2pi,angle2pi>0.0)

        if (angle2pi<halfpi) then
           cosx=cos_32s(angle2pi)
           sinx=sinmultiplier*sqrt(1.0-cosx*cosx)
        elseif(angle2pi<pi) then
           cosx=-cos_32s(pi-angle2pi);
           sinx=sinmultiplier*sqrt(1.0-cosx*cosx)
        elseif (angle2pi<threehalfpi) then
           cosx=-cos_32s(angle2pi-pi);
           sinx=sinmultiplier*sqrt(1.0-cosx*cosx)
        else
           cosx=cos_32s(twopi-angle2pi);
           sinx=sinmultiplier*sqrt(1.0-cosx*cosx)
        endif
    end subroutine sincos32


end module fast_sincos
