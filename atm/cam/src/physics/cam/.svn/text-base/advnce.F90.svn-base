
subroutine advnce(phys_state, cam_out, pbuf)
!-----------------------------------------------------------------------------------
!
! Purpose: The place for parameterizations to call per timestep initializations.
!          Generally this is used to update time interpolated fields from boundary
!          datasets.
!
!-----------------------------------------------------------------------------------
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use chemistry,           only: chem_timestep_init
  use chem_surfvals,       only: chem_surfvals_set
  use ppgrid,              only: begchunk, endchunk
  use physics_types,       only: physics_state
  use phys_buffer,         only: pbuf_size_max, pbuf_fld
  use ghg_data,            only: ghg_data_timestep_init
  use cam3_aero_data,      only: cam3_aero_data_on, cam3_aero_data_timestep_init
  use cam3_ozone_data,     only: cam3_ozone_data_on, cam3_ozone_data_timestep_init
  use radiation,           only: radiation_do
  use tracers,             only: tracers_timestep_init
  use aoa_tracers,         only: aoa_tracers_timestep_init
  use vertical_diffusion,  only: vertical_diffusion_ts_init
  use radheat,             only: radheat_timestep_init
  use solar_data,          only: solar_data_advance
#if ( defined WACCM_PHYS )
  use efield,              only: get_efield
  use iondrag,             only: do_waccm_ions
  use qbo,                 only: qbo_timestep_init
#endif
  use perf_mod

  use prescribed_ozone,    only: prescribed_ozone_adv
  use prescribed_ghg,      only: prescribed_ghg_adv
  use prescribed_aero,     only: prescribed_aero_adv
  use aerodep_flx,         only: aerodep_flx_adv
  use aircraft_emit,       only: aircraft_emit_adv
  use prescribed_volcaero, only: prescribed_volcaero_adv
  use camsrfexch_types,    only: cam_out_t     

  implicit none

  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out
  type(pbuf_fld),      intent(inout), dimension(pbuf_size_max)     :: pbuf  ! physics buffer

  !-----------------------------------------------------------------------------

  ! Chemistry surface values
  call chem_surfvals_set(phys_state)

  ! Solar irradiance
  call solar_data_advance()

  ! Time interpolate for chemistry.
  call chem_timestep_init(phys_state)

  ! Prescribed tracers
  call prescribed_ozone_adv(phys_state)
  call prescribed_ghg_adv(phys_state)
  call prescribed_aero_adv(phys_state)
  call aircraft_emit_adv(phys_state)
  call prescribed_volcaero_adv(phys_state)

  ! prescribed aerosol deposition fluxes
  call aerodep_flx_adv(phys_state, cam_out)

  ! CAM3 prescribed aerosol masses
  if (cam3_aero_data_on) call cam3_aero_data_timestep_init(pbuf, phys_state)

  ! CAM3 prescribed ozone data
  if (cam3_ozone_data_on) call cam3_ozone_data_timestep_init(pbuf, phys_state)

  ! Time interpolate data models of gasses in pbuf
  call ghg_data_timestep_init(pbuf, phys_state)

  ! Upper atmosphere radiative processes
  call radheat_timestep_init(phys_state)
 
  ! Time interpolate for vertical diffusion upper boundary condition
  call vertical_diffusion_ts_init(phys_state)

#if ( defined WACCM_PHYS )
  if (do_waccm_ions) then
     ! Compute the electric field
     call t_startf ('efield')
     call get_efield
     call t_stopf ('efield')
  endif
  !----------------------------------------------------------------------
  ! update QBO data for this time step
  !----------------------------------------------------------------------
  call qbo_timestep_init
#endif

  ! Time interpolate for tracers, if appropriate
  call tracers_timestep_init(phys_state)

  ! age of air tracers
  call aoa_tracers_timestep_init(phys_state)

end subroutine advnce
