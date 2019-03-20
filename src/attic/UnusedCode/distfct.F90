!*****************************************************************************
!                    Distribution Function Routines
!                               M. Shay
!                       Started 3/16/2006
!
!*****************************************************************************

#include "param"



!-----------------------------------------------------------------------------
!            Calculate distribution function
!-----------------------------------------------------------------------------

Subroutine calcdist
#ifdef outputdist
! M. Shay, 3/17/2006. You are asking right now, "Why the hell did Mike rename
! all the variables with a 1 after them?" Why even pass them to
! calcdist as opposed to just using the modules? It turns out that if
! I just used the module and don't pass them, every variable that has
! kind=4 (single precision) is zero unless I surround it with real(,kind=8)
! I don't want to do that, so I have to pass it as an argument.

  Use pe_env
  Use partfield
  Use timer
  Implicit none
  Real(kind=prk), Parameter :: epsilon=1e-4
  Integer :: i,ix,iy,iz,ixpos,lvxpos
  Real(kind=drk), Parameter:: dvxi=2.*vxmaxi/(nvxi-1.) &
                , dvxe=2.*vxmaxe/(nvxe-1.) &
                , dx=lx/(nx*pex)
  Logical :: non_periodic

  disti=0. ; diste=0.

!  call start_timer()

  do i=1,np_i
    ixpos=nint(rvi(1,i)/dx - my_pex*nx + 0.5)
    lvxpos = nint(1. + (rvi(4,i) + vxmaxi)/dvxi)

    if ( (lvxpos >=1) .and. (lvxpos <= nvxi) ) then
      disti(ixpos,lvxpos) = disti(ixpos,lvxpos) + 1.
    end if
  end do

  do i=1,np_e
    ixpos=nint(rve(1,i)/dx - my_pex*nx + 0.5)
    lvxpos = nint(1. + (rve(4,i) + vxmaxe)/dvxe)
    
    if ( (lvxpos >=1) .and. (lvxpos <= nvxe) ) then
      diste(ixpos,lvxpos) = diste(ixpos,lvxpos) + 1.
    end if
  end do
#endif
End Subroutine calcdist


