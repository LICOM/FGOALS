module buffer

!----------------------------------------------------------------------- 
! 
! Purpose: 
!   Definition and initialization of time-cycling physics arrays
!
! Author: 
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst
  use ppgrid, only: pcols, pver, begchunk, endchunk
  use infnan
 
  implicit none
  public
  save

!  real(r8), allocatable :: qrs(:,:,:)      ! shortwave radiative heating rate 
!  real(r8), allocatable :: qrl(:,:,:)      ! longwave  radiative heating rate 

  real(r8), allocatable :: pblht(:,:)      ! PBL height
  real(r8), allocatable :: tpert(:,:)      ! temperature perturbation (PBL)
  real(r8), allocatable :: qpert(:,:,:)    ! moisture/constituent perturb.(PBL)

CONTAINS

  subroutine initialize_buffer
!
! Allocate memory
!
!    allocate (qrs   (pcols,pver,begchunk:endchunk))     
!    allocate (qrl   (pcols,pver,begchunk:endchunk))     
    allocate (pblht (pcols,begchunk:endchunk))       
    allocate (tpert (pcols,begchunk:endchunk))       
    allocate (qpert (pcols,pcnst,begchunk:endchunk)) 
!
! Initialize to NaN or Inf
!
!    qrs   (:,:,:) = inf
!    qrl   (:,:,:) = inf
    pblht (:,:) = inf
    tpert (:,:) = inf
    qpert (:,:,:) = inf

    return
  end subroutine initialize_buffer

end module buffer
