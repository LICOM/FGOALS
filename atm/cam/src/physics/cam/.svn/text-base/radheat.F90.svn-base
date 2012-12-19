
module radheat
!-----------------------------------------------------------------------
!
! Purpose:  Provide an interface to convert shortwave and longwave
!           radiative heating terms into net heating.
!
!           This module provides a hook to allow incorporating additional
!           radiative terms (eUV heating and nonLTE longwave cooling).
! 
! Original version: B.A. Boville
!-----------------------------------------------------------------------

use shr_kind_mod,  only: r8 => shr_kind_r8
use ppgrid,        only: pcols, pver
use physics_types, only: physics_state, physics_ptend, physics_ptend_init
use phys_buffer,   only: pbuf_size_max, pbuf_fld

implicit none
private
save

! Public interfaces
public  &
   radheat_defaultopts,   &!
   radheat_setopts,       &!
   radheat_init,          &!
   radheat_timestep_init, &!
   radheat_tend            ! return net radiative heating

!===============================================================================
contains
!===============================================================================

subroutine radheat_defaultopts( &
   nlte_use_mo_out, itgcmcyc_out, cftgcm_out)

!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   use chemistry, only: chem_is

   logical,          intent(out), optional :: nlte_use_mo_out
   integer,          intent(out), optional :: itgcmcyc_out
   character(len=*), intent(out), optional :: cftgcm_out
!-----------------------------------------------------------------------

end subroutine radheat_defaultopts

!================================================================================================

subroutine radheat_setopts( &
   nlte_use_mo_in, itgcmcyc_in, cftgcm_in)

!----------------------------------------------------------------------- 
! Purpose: Set runtime options
!-----------------------------------------------------------------------

   logical,          intent(in), optional :: nlte_use_mo_in
   integer,          intent(in), optional :: itgcmcyc_in
   character(len=*), intent(in), optional :: cftgcm_in
!-----------------------------------------------------------------------

end subroutine radheat_setopts

!================================================================================================

subroutine radheat_init(hypm)

   use pmgrid, only: plev

   real(r8), intent(in) :: hypm(plev)

end subroutine radheat_init

!================================================================================================

subroutine radheat_timestep_init (state)
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk

    type(physics_state), intent(in):: state(begchunk:endchunk)                 

end subroutine radheat_timestep_init

!================================================================================================

subroutine radheat_tend(state, pbuf, ptend, qrl, qrs, fsns, &
                        fsnt, flns, flnt, asdir, net_flx)
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and the associated net
! boundary flux.
!-----------------------------------------------------------------------

! Arguments
   type(physics_state), intent(in)  :: state             ! Physics state variables
   type(pbuf_fld),      intent(in), dimension(pbuf_size_max) :: pbuf ! only exists for waccm call
   type(physics_ptend), intent(out) :: ptend             ! indivdual parameterization tendencie
   real(r8),            intent(in)  :: qrl(pcols,pver)   ! longwave heating
   real(r8),            intent(in)  :: qrs(pcols,pver)   ! shortwave heating
   real(r8),            intent(in)  :: fsns(pcols)       ! Surface solar absorbed flux
   real(r8),            intent(in)  :: fsnt(pcols)       ! Net column abs solar flux at model top
   real(r8),            intent(in)  :: flns(pcols)       ! Srf longwave cooling (up-down) flux
   real(r8),            intent(in)  :: flnt(pcols)       ! Net outgoing lw flux at model top
   real(r8),            intent(in)  :: asdir(pcols)      ! shortwave, direct albedo
   real(r8),            intent(out) :: net_flx(pcols)  


! Local variables
   integer :: i
   integer :: ncol
!-----------------------------------------------------------------------

   ncol = state%ncol
   call physics_ptend_init(ptend)

   ptend%name       = 'radheat'
   ptend%s(:ncol,:) = (qrs(:ncol,:) + qrl(:ncol,:))
   ptend%ls         = .TRUE.

   do i = 1, ncol
      net_flx(i) = fsnt(i) - fsns(i) - flnt(i) + flns(i)
   end do

end subroutine radheat_tend
end module radheat
