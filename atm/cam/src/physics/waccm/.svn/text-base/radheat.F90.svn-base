
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
  use spmd_utils,    only: masterproc
  use ppgrid,        only: pcols, pver
  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use physconst,     only: gravit, cpair
  use perf_mod
  use cam_logfile,   only: iulog

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

! Namelist variables
  logical :: nlte_use_mo = .true. ! Determines which constituents are used from NLTE calculations
                                  !  = .true. uses prognostic constituents
                                  !  = .false. uses constituents from prescribed dataset waccm_forcing_file

! Private variables for merging heating rates
  integer :: ntop_qrs_cam             ! top level for pure cam solar heating
  integer :: ntop_qrl_cam             ! top level for pure cam long wave heating
  integer :: nbot_qrs_mlt             ! bottom level for pure M/LT solar heating
  integer :: nbot_qrl_mlt             ! bottom level for pure M/LT long wave heating
  real(r8):: qrs_wt(pver)             ! merge weight for cam solar heating
  real(r8):: qrl_wt(pver)             ! merge weight for cam long wave heating

  logical :: waccm_heating

!===============================================================================
contains
!===============================================================================

  subroutine radheat_defaultopts( nlte_use_mo_out )

!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

    use chemistry, only: chem_is

    logical,          intent(out), optional :: nlte_use_mo_out
!-----------------------------------------------------------------------

    if ( present(nlte_use_mo_out) ) then
       nlte_use_mo_out = nlte_use_mo
    end if

  end subroutine radheat_defaultopts

!================================================================================================

  subroutine radheat_setopts( nlte_use_mo_in )

!----------------------------------------------------------------------- 
! Purpose: Set runtime options
!-----------------------------------------------------------------------

    logical,          intent(in), optional :: nlte_use_mo_in
!-----------------------------------------------------------------------

    if ( present(nlte_use_mo_in) ) then
       nlte_use_mo = nlte_use_mo_in
    end if

  end subroutine radheat_setopts

!================================================================================================

  subroutine radheat_init(hypm)

    use pmgrid,           only: plev
    use nlte_lw,          only: nlte_init
    use cam_history,      only: add_default, addfld, phys_decomp

    real(r8), intent(in) :: hypm(plev)

    real(r8) :: psh(pver)                                ! pressure scale height
    integer  :: k
!-----------------------------------------------------------------------

! Initialize merging weights and regions for heating rates
!    ntop_qrs_cam = 0
    ntop_qrs_cam = 1
    nbot_qrs_mlt = 0
!    ntop_qrl_cam = 0
    ntop_qrl_cam = 1
    nbot_qrl_mlt = 0
    qrs_wt(:) = 0._r8
    qrl_wt(:) = 0._r8

! Define pressure scale height (psh)
    do k=1,pver
       psh(k)=log(1e5_r8/hypm(k))
       qrs_wt(k) = 1._r8 - tanh( (psh(k) - 9._r8)/.75_r8 )
       qrl_wt(k) = 1._r8-tanh( (psh(k) - 8.57_r8) / 0.71_r8 )
    end do

    do k=1,pver
       if (psh(k) .gt. 10._r8  ) nbot_qrs_mlt = k
       if (psh(k) .ge.  9._r8  ) ntop_qrs_cam = k+1
       if (psh(k) .gt. 10._r8  ) nbot_qrl_mlt = k
       if (psh(k) .ge.  8.57_r8) ntop_qrl_cam = k+1
    end do
    if (masterproc) then
       write(iulog,*) 'RADHEAT_INIT: NBOT_QRS_MLT, NTOP_QRS_CAM ',  nbot_qrs_mlt, ntop_qrs_cam
       write(iulog,*) 'RADHEAT_INIT: QRS_WT ', qrs_wt(nbot_qrs_mlt+1 : ntop_qrs_cam-1)
       write(iulog,*) 'RADHEAT_INIT: NBOT_QRL_MLT, NTOP_QRL_CAM ',  nbot_qrl_mlt, ntop_qrl_cam
       write(iulog,*) 'RADHEAT_INIT: QRL_WT ', qrl_wt(nbot_qrl_mlt+1 : ntop_qrl_cam-1)
    end if

    waccm_heating = nbot_qrs_mlt > 0 .and. ntop_qrs_cam > 1

    if (waccm_heating) then
       call nlte_init(hypm, nlte_use_mo)
    endif

! Add history variables to master field list
    call addfld ('QRL_TOT  ','K/s     ',plev, 'A','Merged LW heating: QRL+QRLNLTE',phys_decomp)
    call addfld ('QRS_TOT  ','K/s     ',plev, 'A','Merged SW heating: QRS+QCP+QRS_EUV+QRS_CO2NIR+QRS_AUR+QTHERMAL',phys_decomp)

    call addfld ('QRS_TOT_24_COS','K/s  ',plev, 'A','SW heating 24hr. cos coeff.',phys_decomp)
    call addfld ('QRS_TOT_24_SIN','K/s  ',plev, 'A','SW heating 24hr. sin coeff.',phys_decomp)
    call addfld ('QRS_TOT_12_COS','K/s  ',plev, 'A','SW heating 12hr. cos coeff.',phys_decomp)
    call addfld ('QRS_TOT_12_SIN','K/s  ',plev, 'A','SW heating 12hr. sin coeff.',phys_decomp)

! Add default history variables to files
    call add_default ('QRL_TOT', 1, ' ')
    call add_default ('QRS_TOT', 1, ' ')
    call add_default ('QRS_TOT', 2, ' ')

  end subroutine radheat_init

!================================================================================================

  subroutine radheat_timestep_init (state)

    use nlte_lw, only: nlte_timestep_init
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk

    type(physics_state), intent(in):: state(begchunk:endchunk)                 

    if (waccm_heating) then
       call nlte_timestep_init (state)
    endif

  end subroutine radheat_timestep_init

!================================================================================================

  subroutine radheat_tend(state, pbuf, ptend, qrl, qrs, fsns, &
       fsnt, flns, flnt, asdir, net_flx)
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and the associated net
! boundary flux.
!
! This routine provides the waccm hook for computing nonLTE cooling and 
! eUV heating. 
!-----------------------------------------------------------------------

    use cam_history,       only: outfld
    use nlte_lw,           only: nlte_tend
    use mo_waccm_hrates,   only: waccm_hrates, has_hrates
    use waccm_forcing,     only: get_solar
    use phys_buffer,       only: pbuf_size_max, pbuf_fld
    use tidal_diag,        only: get_tidal_coeffs

! Arguments
    type(physics_state), intent(in)  :: state             ! Physics state variables
    type(pbuf_fld), intent(in), dimension(pbuf_size_max) :: pbuf ! Source for Trop. O3?
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
    integer  :: i, k
    integer  :: ncol                                ! number of atmospheric columns
    integer :: lchnk                                ! chunk identifier
    real(r8) :: qrl_mrg(pcols,pver)                 ! merged LW heating
    real(r8) :: qrl_mlt(pcols,pver)                 ! M/LT longwave heating rates
    real(r8) :: qrs_mrg(pcols,pver)                 ! merged SW heating
    real(r8) :: qrs_mlt(pcols,pver)                 ! M/LT solar heating rates
    real(r8) :: qout(pcols,pver)                    ! temp for outfld call
    real(r8) :: dcoef(4)                            ! for tidal component of heating
!-----------------------------------------------------------------------

    ncol  = state%ncol
    lchnk = state%lchnk
    call physics_ptend_init(ptend)

! WACCM interactive heating rate
    if (waccm_heating) then
       call t_startf( 'hrates' )
       if (has_hrates) then
          call waccm_hrates(ncol, state, asdir, qrs_mlt)
       else
          call get_solar(ncol, lchnk, qrs_mlt)
       endif
       call t_stopf( 'hrates' )
    else
       qrs_mlt(:,:) = 0._r8
    endif

! Merge cam solar heating for lower atmosphere with M/LT heating
    call merge_qrs (ncol, qrs, qrs_mlt, qrs_mrg)
    qout(:ncol,:) = qrs_mrg(:ncol,:)/cpair
    call outfld ('QRS_TOT', qout, pcols, lchnk)

! Output tidal coefficients of total SW heating
    call get_tidal_coeffs( dcoef )
    call outfld( 'QRS_TOT_24_SIN', qout(:ncol,:)*dcoef(1), ncol, lchnk )
    call outfld( 'QRS_TOT_24_COS', qout(:ncol,:)*dcoef(2), ncol, lchnk )
    call outfld( 'QRS_TOT_12_SIN', qout(:ncol,:)*dcoef(3), ncol, lchnk )
    call outfld( 'QRS_TOT_12_COS', qout(:ncol,:)*dcoef(4), ncol, lchnk )

    if (waccm_heating) then
       call t_startf( 'nltedrv' )
       call nlte_tend(state, pbuf, qrl_mlt)
       call t_stopf( 'nltedrv' )
    else
       qrl_mlt(:,:) = 0._r8
    endif

!   Merge cam long wave heating for lower atmosphere with M/LT (nlte) heating
    call merge_qrl (ncol, qrl, qrl_mlt, qrl_mrg)
    qout(:ncol,:) = qrl_mrg(:ncol,:)/cpair
    call outfld ('QRL_TOT', qout, pcols, lchnk)

    ptend%name       = 'radheat'
    ptend%s(:ncol,:) = qrs_mrg(:ncol,:) + qrl_mrg(:ncol,:)
    ptend%ls         = .TRUE.

    net_flx = 0._r8
    do k = 1, pver
       do i = 1, ncol
          net_flx(i) = net_flx(i) + ptend%s(i,k)*state%pdel(i,k)/gravit
       end do
    end do

  end subroutine radheat_tend

!================================================================================================

  subroutine merge_qrs (ncol, hcam, hmlt, hmrg)
!
!  Merges short wave heating rates
!
    implicit none

!-----------------Input arguments-----------------------------------
    integer ncol

    real(r8), intent(in)  :: hmlt(pcols,pver)                ! Upper atmosphere heating rates
    real(r8), intent(in)  :: hcam(pcols,pver)                ! CAM heating rate
    real(r8), intent(out) :: hmrg(pcols,pver)                ! merged heating rates

!-----------------Local workspace------------------------------------

    integer k

!--------------------------------------------------------------------

! Pure M/LT region
    do k = 1, nbot_qrs_mlt
       hmrg(:ncol,k) = cpair*hmlt(:ncol,k)
    end do

! Merge region
    do k = nbot_qrs_mlt+1, ntop_qrs_cam-1
       hmrg(:ncol,k) = qrs_wt(k)*hcam(:ncol,k) + (1._r8 - qrs_wt(k))*cpair*hmlt(:ncol,k) 
    end do

! Pure cam region
    do k = ntop_qrs_cam, pver
       hmrg(:ncol,k) = hcam(:ncol,k)    
    end do

  end subroutine merge_qrs

!==================================================================================================

  subroutine merge_qrl (ncol, hcam, hmlt, hmrg)
!
!  Merges long wave heating rates
!
!-----------------Input arguments-----------------------------------
    integer ncol

    real(r8), intent(in)  :: hmlt(pcols,pver)               ! Upper atmosphere heating rates
    real(r8), intent(in)  :: hcam(pcols,pver)               ! CAM heating rate
    real(r8), intent(out) :: hmrg(pcols,pver)               ! merged heating rates

!-----------------Local workspace------------------------------------

    integer k

!--------------------------------------------------------------------

! Pure M/LT region
    do k = 1, nbot_qrl_mlt
       hmrg(:ncol,k) = hmlt(:ncol,k)
    end do

! Merge region
    do k = nbot_qrl_mlt+1, ntop_qrl_cam-1
       hmrg(:ncol,k) = qrl_wt(k) * hcam(:ncol,k) + (1._r8-qrl_wt(k)) * hmlt(:ncol,k) 
    end do

! Pure cam region
    do k = ntop_qrl_cam, pver
       hmrg(:ncol,k)   = hcam(:ncol,k)    
    end do

  end subroutine merge_qrl

end module radheat
