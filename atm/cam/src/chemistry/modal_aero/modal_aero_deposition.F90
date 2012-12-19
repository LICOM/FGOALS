module modal_aero_deposition

!------------------------------------------------------------------------------------------------
! Purpose:
!
! Partition the contributions from modal components of wet and dry 
! deposition at the surface into the fields passed to the coupler.
!
! *** N.B. *** Currently only a simple scheme for the 3-mode version
!              of MAM has been implemented.
!
! Revision history:
! Feb 2009  M. Flanner, B. Eaton   Original version for trop_mam3.
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use camsrfexch_types, only: cam_out_t     
use constituents,     only: pcnst, cnst_get_ind
use ppgrid,           only: pcols

implicit none
private
save

public :: &
   modal_aero_deposition_init, &
   set_srf_drydep,             &
   set_srf_wetdep

! Private module data
integer :: idx_bc1  = -1
integer :: idx_pom1 = -1
integer :: idx_soa1 = -1
integer :: idx_soa2 = -1
integer :: idx_dst1 = -1
integer :: idx_dst3 = -1
integer :: idx_ncl3 = -1
integer :: idx_so43 = -1
integer :: idx_num3 = -1

!==============================================================================
contains
!==============================================================================

subroutine modal_aero_deposition_init()

! set aerosol indices for re-mapping surface deposition fluxes:
! *_a1 = accumulation mode
! *_a2 = aitken mode
! *_a3 = coarse mode

   ! Currently only trop_mam3 scheme is implemented.
#ifndef MODAL_AERO_3MODE  
   return
#endif

   call cnst_get_ind('bc_a1',  idx_bc1)
   call cnst_get_ind('pom_a1', idx_pom1)
   call cnst_get_ind('soa_a1', idx_soa1)
   call cnst_get_ind('soa_a2', idx_soa2)
   call cnst_get_ind('dst_a1', idx_dst1)
   call cnst_get_ind('dst_a3', idx_dst3)
   call cnst_get_ind('ncl_a3', idx_ncl3)
   call cnst_get_ind('so4_a3', idx_so43)
   call cnst_get_ind('num_a3', idx_num3)

end subroutine modal_aero_deposition_init

!==============================================================================
subroutine set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)

! Set surface wet deposition fluxes passed to coupler.

   ! Arguments:
   real(r8), intent(in) :: aerdepwetis(pcols,pcnst)  ! aerosol wet deposition (interstitial)
   real(r8), intent(in) :: aerdepwetcw(pcols,pcnst)  ! aerosol wet deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i
   integer :: ncol                      ! number of columns
   !----------------------------------------------------------------------------

   ! Currently only trop_mam3 scheme is implemented.
#ifndef MODAL_AERO_3MODE  
   return
#endif

   ncol = cam_out%ncol

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface, 
   !        dry deposition fluxes are positive into surface.
   !        CLM wants positive definite fluxes.
   do i = 1, ncol
      ! black carbon fluxes
      cam_out%bcphiwet(i) = -(aerdepwetis(i,idx_bc1)+aerdepwetcw(i,idx_bc1))

      ! organic carbon fluxes
      cam_out%ocphiwet(i) = -(aerdepwetis(i,idx_pom1)+aerdepwetis(i,idx_soa1)+aerdepwetcw(i,idx_pom1)+aerdepwetcw(i,idx_soa1))

      ! dust fluxes
      !
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      cam_out%dstwet1(i) = -(aerdepwetis(i,idx_dst1)+aerdepwetcw(i,idx_dst1))
      
      !  A. Simple: Assign all coarse-mode dust to bulk size bin 3:
      cam_out%dstwet2(i) = 0._r8
      cam_out%dstwet3(i) = -(aerdepwetis(i,idx_dst3)+aerdepwetcw(i,idx_dst3))
      cam_out%dstwet4(i) = 0._r8

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphiwet(i) .lt. 0._r8) cam_out%bcphiwet(i) = 0._r8
      if (cam_out%ocphiwet(i) .lt. 0._r8) cam_out%ocphiwet(i) = 0._r8
      if (cam_out%dstwet1(i)  .lt. 0._r8) cam_out%dstwet1(i)  = 0._r8
      if (cam_out%dstwet3(i)  .lt. 0._r8) cam_out%dstwet3(i)  = 0._r8
   enddo

end subroutine set_srf_wetdep

!==============================================================================

subroutine set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)

! Set surface dry deposition fluxes passed to coupler.
   
   ! Arguments:
   real(r8), intent(in) :: aerdepdryis(pcols,pcnst)  ! aerosol dry deposition (interstitial)
   real(r8), intent(in) :: aerdepdrycw(pcols,pcnst)  ! aerosol dry deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i
   integer :: ncol                      ! number of columns
   !----------------------------------------------------------------------------

   ! Currently only trop_mam3 scheme is implemented.
#ifndef MODAL_AERO_3MODE  
   return
#endif

   ncol = cam_out%ncol

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface, 
   !        dry deposition fluxes are positive into surface.
   !        CLM wants positive definite fluxes.
   do i = 1, ncol
      ! black carbon fluxes
      cam_out%bcphidry(i) = aerdepdryis(i,idx_bc1)+aerdepdrycw(i,idx_bc1)
      cam_out%bcphodry(i) = 0._r8

      ! organic carbon fluxes
      cam_out%ocphidry(i) = aerdepdryis(i,idx_pom1)+aerdepdryis(i,idx_soa1)+aerdepdrycw(i,idx_pom1)+aerdepdrycw(i,idx_soa1)
      cam_out%ocphodry(i) = aerdepdryis(i,idx_soa2)+aerdepdrycw(i,idx_soa2)

      ! dust fluxes
      !
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      cam_out%dstdry1(i) = aerdepdryis(i,idx_dst1)+aerdepdrycw(i,idx_dst1)
      
      ! Two options for partitioning deposition into bins 2-4:
      !  A. Simple: Assign all coarse-mode dust to bulk size bin 3:
      cam_out%dstdry2(i) = 0._r8
      cam_out%dstdry3(i) = aerdepdryis(i,idx_dst3)+aerdepdrycw(i,idx_dst3)
      cam_out%dstdry4(i) = 0._r8

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphidry(i) .lt. 0._r8) cam_out%bcphidry(i) = 0._r8
      if (cam_out%bcphodry(i) .lt. 0._r8) cam_out%bcphodry(i) = 0._r8
      if (cam_out%ocphidry(i) .lt. 0._r8) cam_out%ocphidry(i) = 0._r8
      if (cam_out%ocphodry(i) .lt. 0._r8) cam_out%ocphodry(i) = 0._r8
      if (cam_out%dstdry1(i)  .lt. 0._r8) cam_out%dstdry1(i)  = 0._r8
      if (cam_out%dstdry3(i)  .lt. 0._r8) cam_out%dstdry3(i)  = 0._r8
   enddo

end subroutine set_srf_drydep


!==============================================================================

end module modal_aero_deposition
