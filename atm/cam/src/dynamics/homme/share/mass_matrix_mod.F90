#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module mass_matrix_mod
  use kinds, only : real_kind
  use dimensions_mod, only : np, nelemd, nv
  use quadrature_mod, only : quadrature_t, gauss ,gausslobatto
  use element_mod, only : element_t
  use parallel_mod, only : parallel_t
  use edge_mod, only : edgebuffer_t,edgevpack,edgevunpack, &
       freeedgebuffer,initedgebuffer  
  use bndry_mod, only : bndry_exchangev
implicit none
private

  public :: mass_matrix

contains

! ===========================================
! mass_matrix:
!
! Compute the mass matrix for each element...
! ===========================================

  subroutine mass_matrix(par,elem)

    type (parallel_t),intent(in) :: par
    type (element_t) :: elem(:)

    type (EdgeBuffer_t)    :: edge

    real(kind=real_kind)  da                     ! area element

    type (quadrature_t) :: gv
    type (quadrature_t) :: gp

    integer ii
    integer i,j
    integer kptr
    integer iptr

    ! ===================
    ! begin code
    ! ===================

    call initEdgeBuffer(edge,1)

    ! =================================================
    ! mass matrix on the velocity grid
    ! =================================================    

    gv=gausslobatto(nv)
 
    do ii=1,nelemd
       da=0.25D0*elem(ii)%dx * elem(ii)%dy
       do j=1,nv
          do i=1,nv
             elem(ii)%mv(i,j)=da*(gv%weights(i)*gv%weights(j))
             elem(ii)%rmv(i,j)=elem(ii)%mv(i,j)
          end do
       end do

       kptr=0
       call edgeVpack(edge,elem(ii)%rmv,1,kptr,elem(ii)%desc)

    end do

    ! ==============================
    ! Insert boundary exchange here
    ! ==============================

    call bndry_exchangeV(par,edge)

    do ii=1,nelemd

       kptr=0
       call edgeVunpack(edge,elem(ii)%rmv,1,kptr,elem(ii)%desc)

       do j=1,nv
          do i=1,nv
             elem(ii)%rmv(i,j)=1.0D0/elem(ii)%rmv(i,j)
          end do
       end do

    end do
!$OMP BARRIER

    deallocate(gv%points)
    deallocate(gv%weights)

    ! =============================================
    ! compute spherical element mass matrix
    ! =============================================
    do ii=1,nelemd
       do j=1,nv
          do i=1,nv
             elem(ii)%spheremv(i,j)=elem(ii)%mv(i,j)*elem(ii)%metdet(i,j)
             elem(ii)%rspheremv(i,j)=elem(ii)%spheremv(i,j)
          end do
       end do
       kptr=0
       call edgeVpack(edge,elem(ii)%rspheremv,1,kptr,elem(ii)%desc)
    end do
    call bndry_exchangeV(par,edge)
    do ii=1,nelemd
       kptr=0
       call edgeVunpack(edge,elem(ii)%rspheremv,1,kptr,elem(ii)%desc)
       do j=1,nv
          do i=1,nv
             elem(ii)%rspheremv(i,j)=1.0D0/elem(ii)%rspheremv(i,j)
          end do
       end do
    end do
!$OMP BARRIER



    ! =============================================
    ! compute the mass matrix on the pressure grid
    ! =============================================

    if (np == nv) then
       do ii=1,nelemd
          iptr=1
          do j=1,np
             do i=1,np
                elem(ii)%mp(i,j)=elem(ii)%mv(i,j)
                iptr=iptr+1
             end do
          end do
       end do
    else if (np == nv-2) then
       gp=gauss(np)
       do ii=1,nelemd
          da=0.25D0*elem(ii)%dx * elem(ii)%dy
          iptr=1
          do j=1,np
             do i=1,np
                elem(ii)%mp(i,j)=da*(gp%weights(i)*gp%weights(j))
                iptr=iptr+1
             end do
          end do
       end do
       deallocate(gp%points)
       deallocate(gp%weights)
    end if
   
    call FreeEdgeBuffer(edge)
       
  end subroutine mass_matrix

end module mass_matrix_mod

