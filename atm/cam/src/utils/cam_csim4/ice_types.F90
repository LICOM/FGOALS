module ice_types

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use ppgrid,        only: pcols, begchunk, endchunk
  implicit none
  public

  type frac_t
     real(r8) :: land(pcols)
  end type frac_t

  type ice_in_t
     real(r8) :: tsocn(pcols)        ! ocean layer temperature
     real(r8) :: tbot(pcols)         ! bot level temperature
     real(r8) :: zbot(pcols)         ! bot level height above surface
     real(r8) :: ubot(pcols)         ! bot level u wind
     real(r8) :: vbot(pcols)         ! bot level v wind
     real(r8) :: qbot(pcols)         ! bot level specific humidity
     real(r8) :: pbot(pcols)         ! bot level pressure
     real(r8) :: flwds(pcols)        !
     real(r8) :: snow(pcols)         ! 
     real(r8) :: soll(pcols)         ! 
     real(r8) :: sols(pcols)         ! 
     real(r8) :: solld(pcols)        !
     real(r8) :: solsd(pcols)        !
     real(r8) :: thbot(pcols)        ! 
     real(r8) :: frzmlt(pcols)       ! 
  end type ice_in_t

  type ice_out_t
     real(r8) :: areafrac(pcols)     ! area fraction of ice over all gridcell
     real(r8) :: aice(pcols)         ! area fraction of ice over open ocean
     real(r8) :: sicthk(pcols)       ! sea ice thickness (m)            
     real(r8) :: asdir(pcols)        ! albedo: shortwave, direct
     real(r8) :: asdif(pcols)        ! albedo: shortwave, diffuse
     real(r8) :: aldir(pcols)        ! albedo: longwave, direct
     real(r8) :: aldif(pcols)        ! albedo: longwave, diffuse
     real(r8) :: lwup(pcols)         ! longwave up radiative flux
     real(r8) :: lhf(pcols)          ! latent heat flux
     real(r8) :: shf(pcols)          ! sensible heat flux
     real(r8) :: wsx(pcols)          ! surface u-stress (N)
     real(r8) :: wsy(pcols)          ! surface v-stress (N)
     real(r8) :: tref(pcols)         ! ref height surface air temp
     real(r8) :: ts(pcols)           ! sfc temp (merged w/ocean if coupled)
     real(r8) :: cflx(pcols)         ! constituent flux (evap)
     real(r8) :: focn(pcols)         ! ocean ice heat flux for basl and lateral melt (<0)
     real(r8) :: fswabs(pcols)       ! SW absorbed in ice
  end type ice_out_t

contains

  subroutine ice_types_alloc(ice_in, ice_out)

    type(ice_in_t) , pointer :: ice_in(:)
    type(ice_out_t), pointer :: ice_out(:)

    integer :: c	

    allocate (ice_in(begchunk:endchunk))

    do c = begchunk,endchunk
       ice_in(c)%tsocn(:)   = 0._r8
       ice_in(c)%tbot(:)    = 0._r8
       ice_in(c)%zbot(:)    = 0._r8
       ice_in(c)%ubot(:)    = 0._r8
       ice_in(c)%vbot(:)    = 0._r8
       ice_in(c)%qbot(:)    = 0._r8
       ice_in(c)%pbot(:)    = 0._r8
       ice_in(c)%flwds(:)   = 0._r8 
       ice_in(c)%snow(:)    = 0._r8
       ice_in(c)%soll(:)    = 0._r8      
       ice_in(c)%sols(:)    = 0._r8      
       ice_in(c)%solld(:)   = 0._r8     
       ice_in(c)%solsd(:)   = 0._r8     
       ice_in(c)%thbot(:)   = 0._r8
       ice_in(c)%frzmlt(:)  = 0._r8
    end do

    allocate (ice_out(begchunk:endchunk))

    do c = begchunk,endchunk
       ice_out(c)%areafrac(:) = 0.0_r8
       ice_out(c)%aice(:)     = 0.0_r8
       ice_out(c)%sicthk(:)   = 0.0_r8
       ice_out(c)%asdir(:)    = 0.0_r8
       ice_out(c)%asdif(:)    = 0.0_r8
       ice_out(c)%aldir(:)    = 0.0_r8
       ice_out(c)%aldif(:)    = 0.0_r8
       ice_out(c)%lwup(:)     = 0.0_r8
       ice_out(c)%lhf(:)      = 0.0_r8
       ice_out(c)%shf(:)      = 0.0_r8
       ice_out(c)%cflx(:)     = 0.0_r8
       ice_out(c)%wsx(:)      = 0.0_r8
       ice_out(c)%wsy(:)      = 0.0_r8
       ice_out(c)%tref(:)     = 0.0_r8
       ice_out(c)%ts(:)       = 0.0_r8
       ice_out(c)%focn(:)     = 0.0_r8
       ice_out(c)%fswabs(:)   = 0.0_r8
    end do

  end subroutine ice_types_alloc

end module ice_types
