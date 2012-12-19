
module aerosol_intr

!---------------------------------------------------------------------------------
! Module to interface the aerosol parameterizations with CAM
! Phil Rasch, Jan 2003
!---------------------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
use spmd_utils,   only: masterproc
use camsrfexch_types, only: cam_in_t, cam_out_t
use ppgrid,       only: pcols, pver, pverp
use physconst,    only: mwdry, mwh2o, gravit
use phys_control, only: cam_physpkg_is
use constituents, only: pcnst, cnst_name
use abortutils,   only: endrun
use perf_mod
use cam_logfile,  only: iulog

implicit none
private
save

! Public interfaces

public ::&
   aerosol_readnl,                &  ! read aerosol_nl namelist group
   aerosol_register_cnst,         &  ! register consituents
   aerosol_implements_cnst,       &  ! returns true if consituent is implemented by this package
   aerosol_init_cnst,             &  ! initialize mixing ratios if not read from initial file
   prognostic_aerosol_initialize, &  ! initialize (history) variables
   aerosol_drydep_intr,           &  ! interface to dry deposition
   aerosol_wet_intr,              &  ! interface to wet deposition
   aerosol_emis_intr                 ! interface to surface emissions

! Set this flag to .TRUE. to turn on prognostic sea salt
logical, parameter :: def_progsslt = .FALSE.  ! default
logical            :: progsslt = def_progsslt

! dust
! Set this flag to .TRUE. to turn on dust
logical, parameter :: def_dust = .FALSE.  ! default
logical            :: dust = def_dust

! It is useful to know if any of the aerosols are running
! (for instance, to initialize dry deposition module
logical :: is_any_aerosol = .false.

! Physics buffer indices
integer  ::   dgnum_idx            = 0
integer  ::   dgnumwet_idx         = 0
integer  ::   wetdens_ap_idx       = 0
integer  ::   qaerwat_idx          = 0
integer  ::   fracis_idx           = 0

integer  ::    cld_idx             = 0
integer  ::    qme_idx             = 0 
integer  ::    prain_idx           = 0 
integer  ::    nevapr_idx          = 0 
integer  ::    icwmrdp_idx         = 0 
integer  ::    rprddp_idx          = 0 
integer  ::    icwmrsh_idx         = 0 
integer  ::    rprdsh_idx          = 0 
integer  ::    sh_frac_idx         = 0 
integer  ::    dp_frac_idx         = 0 
integer  ::    nevapr_shcu_idx     = 0 
integer  ::    nevapr_dpcu_idx     = 0 

integer  ::    rate1_cw2pr_st_idx  = 0  


! Namelist variables
real(r8)      :: dust_emis_fact = -1.e36   ! tuning parameter for dust emissions
character(cl) :: soil_erod = 'soil_erod'   ! full pathname for soil erodibility dataset


!===============================================================================
contains
!===============================================================================

subroutine aerosol_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'aerosol_readnl'

   namelist /aerosol_nl/ dust_emis_fact, soil_erod

   !-----------------------------------------------------------------------------

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'aerosol_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, aerosol_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(dust_emis_fact,               1, mpir8,   0, mpicom)
   call mpibcast(soil_erod,       len(soil_erod), mpichar, 0, mpicom)
#endif

end subroutine aerosol_readnl

!===============================================================================

  subroutine aerosol_register_cnst
!----------------------------------------------------------------------- 
! 
! Purpose: register aerosols
! 
! Method: 

! Author: P. J. Rasch
! 
!-----------------------------------------------------------------------
    use dust_intr,        only: dust_register_cnst
    use dust_intr,        only: ndst=>dust_number, dust_names
    use progseasalts_intr,only: progseasalts_register_cnst
    use progseasalts_intr,only: nsst, progseasalts_names
    use phys_buffer,      only: pbuf_times, pbuf_add
    use camsrfexch_types, only: hub2atm_setopts
#if ( defined MODAL_AERO )
    use modal_aero_data,  only: ntot_amode
#endif

    use mo_chem_utls, only: get_spc_ndx

    implicit none
!---------------------------Local workspace-----------------------------
    integer :: m                               ! tracer index

    integer :: dst_idx(ndst), slt_idx(nsst)

!-----------------------------------------------------------------------

    do m=1,ndst
       dst_idx(m) = get_spc_ndx( dust_names(m) )
    enddo
    do m=1,nsst
       slt_idx(m) = get_spc_ndx( progseasalts_names(m) )
    enddo

    dust      = all( dst_idx > 0 )
    progsslt  = all( slt_idx > 0 )

    if ( dust ) then
       is_any_aerosol = .true.
    endif

    if ( progsslt ) then
       is_any_aerosol = .true.
    endif
    if ( dust .or. progsslt ) then
       ! tell camsrfexch_types to allocate fv & ram1 -- needed by prodsslts and dust
       call hub2atm_setopts(aero_dust_in=.true.)
    endif

#if ( defined MODAL_AERO ) 
    is_any_aerosol = .true.
    call hub2atm_setopts(aero_dust_in=.true.)
    call pbuf_add( 'DGNUM',      'global', 1, pver, ntot_amode,  dgnum_idx )    
    call pbuf_add( 'DGNUMWET',   'global', 1, pver, ntot_amode,  dgnumwet_idx )
    call pbuf_add( 'WETDENS_AP', 'physpkg', 1, pver, ntot_amode, wetdens_ap_idx )
    call pbuf_add( 'QAERWAT',    'physpkg', 1, pver, ntot_amode, qaerwat_idx )
#endif 

! Request physics buffer space for fields that don't persist across timesteps.
    call pbuf_add('FRACIS' , 'physpkg', 1,pver, pcnst, fracis_idx)

    return
  end subroutine aerosol_register_cnst


!=======================================================================
  function aerosol_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this 
!          package
! 
! Author: T. Henderson
! 
!-----------------------------------------------------------------------
    use dust_intr,         only: dust_implements_cnst
    use progseasalts_intr, only: progseasalts_implements_cnst

     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: aerosol_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     aerosol_implements_cnst = .false.

     if ( dust ) then
        aerosol_implements_cnst = &
             (aerosol_implements_cnst.OR.dust_implements_cnst (name))
     endif
     if ( progsslt ) then
        aerosol_implements_cnst = &
             (aerosol_implements_cnst.OR.progseasalts_implements_cnst (name))
     endif


  end function aerosol_implements_cnst


!=======================================================================
  subroutine aerosol_init_cnst(name, q, gcid)
!----------------------------------------------------------------------- 
! 
! Purpose:  
! Set initial mass mixing ratios.  
!
!-----------------------------------------------------------------------
    use drydep_mod,       only: inidrydep
    use dust_intr,        only: dust_implements_cnst, dust_init_cnst
    use progseasalts_intr,only: progseasalts_implements_cnst, progseasalts_init_cnst


    implicit none
!-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name         ! constituent name
    
    real(r8), intent(out) :: q(:,:)            !  mass mixing ratio
    integer, intent(in) :: gcid(:)             ! global column id
!-----------------------------------------------------------------------

    if ( dust ) then
       if (dust_implements_cnst(name)) then
          call dust_init_cnst(name, q, gcid)
       endif
    endif

    if ( progsslt ) then
       if (progseasalts_implements_cnst(name)) then
          call progseasalts_init_cnst(name, q, gcid)
       endif
    endif

  end subroutine aerosol_init_cnst


!===============================================================================
  subroutine prognostic_aerosol_initialize(phys_state)
!----------------------------------------------------------------------- 
! 
! Purpose: initialize aerosol parameterizations
!          (declare history variables)
! 
! Method: 
! call subroutines
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use physics_types,     only: physics_state
    use dust_intr,         only: dust_initialize
    use progseasalts_intr, only: progseasalts_initialize
    use drydep_mod,        only: inidrydep
    use physconst,         only: rair
    use phys_buffer,       only: pbuf_get_fld_idx
    implicit none
    type(physics_state), intent(in) :: phys_state(:)
    
#if ( defined MODAL_AERO )
    call dust_initialize(dust_emis_fact, soil_erod)
    call progseasalts_initialize
#else
    if ( dust ) then
       call dust_initialize(dust_emis_fact, soil_erod)
    endif
    if ( progsslt ) then
       call progseasalts_initialize
    endif
#endif

! Dry deposition needs to be initialized if any of the aerosols
! are running.
    if ( is_any_aerosol ) then
       call inidrydep (rair, gravit )
    endif

    
   cld_idx             = pbuf_get_fld_idx('CLD')    
   qme_idx             = pbuf_get_fld_idx('QME')    
   prain_idx           = pbuf_get_fld_idx('PRAIN')  
   nevapr_idx          = pbuf_get_fld_idx('NEVAPR') 
   icwmrdp_idx         = pbuf_get_fld_idx('ICWMRDP') 
   rprddp_idx          = pbuf_get_fld_idx('RPRDDP')  
   icwmrsh_idx         = pbuf_get_fld_idx('ICWMRSH') 
   rprdsh_idx          = pbuf_get_fld_idx('RPRDSH')  
   sh_frac_idx         = pbuf_get_fld_idx('SH_FRAC' )
   dp_frac_idx         = pbuf_get_fld_idx('DP_FRAC') 
   nevapr_shcu_idx     = pbuf_get_fld_idx('NEVAPR_SHCU') 
   nevapr_dpcu_idx     = pbuf_get_fld_idx('NEVAPR_DPCU') 

#if ( defined MODAL_AERO )
   rate1_cw2pr_st_idx  = pbuf_get_fld_idx( 'RATE1_CW2PR_ST' ) 
#endif
    return
  end subroutine prognostic_aerosol_initialize

!===============================================================================
  subroutine aerosol_wet_intr (state, ptend, dt, pbuf, cam_out, dlf)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to wet processing of all aerosols
! 
! Method: 
!  use a modified version of the scavenging parameterization described in
!     Barth et al, 2000, JGR (sulfur cycle paper)
!     Rasch et al, 2001, JGR (INDOEX paper)
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use cam_history,   only: outfld
    use physics_types, only: physics_state, physics_ptend
    use phys_buffer,   only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx

    use phys_grid,     only: get_lat_all_p, get_lon_all_p, get_rlat_all_p, get_rlon_all_p
    use dust_intr,        only: dust_wet_intr
    use progseasalts_intr,only: progseasalts_wet_intr

    use wetdep,        only: clddiag
    use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday
    use time_manager,  only: is_perpetual, get_nstep
    use scyc,          only: scyc_idx1
    use mz_aerosols_intr,       only:  mz_aero_wet_intr
#if ( defined MODAL_AERO )
    use modal_aero_calcsize,    only:  modal_aero_calcsize_sub
    use modal_aero_wateruptake, only:  modal_aero_wateruptake_sub
    use modal_aero_data,        only:  ntot_amode
#endif
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables

    type(physics_ptend), intent(inout) :: ptend        ! indivdual parameterization tendencies
    type(pbuf_fld),      intent(inout) :: pbuf(pbuf_size_max)
    type(cam_out_t),     intent(inout) :: cam_out      ! export state

    real(r8), intent(in) :: dlf(pcols,pver)            ! shallow+deep convective detrainment [kg/kg/s]


! local vars
    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    real(r8):: rainmr(pcols,pver)  ! mixing ratio of rain within cloud volume
    real(r8):: cldv(pcols,pver)     ! cloudy volume undergoing wet chem and scavenging
    real(r8):: cldvcu(pcols,pver)   ! Convective precipitation area at the top interface of current layer
    real(r8):: cldvst(pcols,pver)   ! Stratiform precipitation area at the top interface of current layer 
    integer ix
    integer  lat(pcols)                  ! latitude indices
    real(r8) clat(pcols)                    ! latitudes
    integer  lon(pcols)                  ! longtitude indices
    real(r8) clon(pcols)                    ! longitudes
    integer  nstep
    real(r8) conicw(pcols,pver)          ! convective in-cloud water
    real(r8) cmfdqr(pcols,pver)          ! convective production of rain
    real(r8) cldc(pcols,pver)            ! convective cloud fraction, currently empty
    real(r8) clds(pcols,pver)            ! Stratiform cloud fraction
    real(r8) evapc(pcols,pver)           ! Evaporation rate of convective precipitation

#if ( defined MODAL_AERO )
    integer  i, k
#endif


! physics buffer 
    integer itim
    real(r8), pointer, dimension(:,:) :: cldn       ! cloud fraction
    real(r8), pointer, dimension(:,:) :: cme
    real(r8), pointer, dimension(:,:) :: prain
    real(r8), pointer, dimension(:,:) :: evapr
    real(r8), pointer, dimension(:,:) ::  icwmrdp ! in cloud water mixing ratio, deep convection
    real(r8), pointer, dimension(:,:) ::  rprddp ! rain production, deep convection
    real(r8), pointer, dimension(:,:) ::  icwmrsh ! in cloud water mixing ratio, deep convection
    real(r8), pointer, dimension(:,:) ::  rprdsh ! rain production, deep convection
    real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
! Dec.29.2009. Sungsu
    real(r8), pointer, dimension(:,:) ::  sh_frac  ! Shallow convective cloud fraction
    real(r8), pointer, dimension(:,:) ::  dp_frac  ! Deep convective cloud fraction
    real(r8), pointer, dimension(:,:) ::  evapcsh  ! Evaporation rate of shallow convective precipitation >=0.
    real(r8), pointer, dimension(:,:) ::  evapcdp  ! Evaporation rate of deep    convective precipitation >=0.
! Dec.29.2009. Sungsu 
#if ( defined MODAL_AERO )
    real(r8), pointer, dimension(:,:,:) :: dgnum_pbuf, dgnumwet_pbuf, wetdens_pbuf
    real(r8), pointer, dimension(:,:,:) :: qaerwat  ! aerosol water
    real(r8), pointer, dimension(:,:) :: rate1ord_cw2pr_st   ! 1st order rate for direct conversion of
                                                    ! strat. cloud water to precip (1/s)    ! rce 2010/05/01
#endif


    integer ncol,lchnk

    nstep = get_nstep()

    calday = get_curr_calday()
    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if
    ncdate = yr*10000 + mon*100 + day

    ncol = state%ncol
    lchnk = state%lchnk

   call get_lat_all_p(lchnk, ncol, lat)
   call get_lon_all_p(lchnk, ncol, lon)
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)

! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   cldn => pbuf(cld_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   cme  => pbuf(qme_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   prain  => pbuf(prain_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   evapr  => pbuf(nevapr_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

   fracis  => pbuf(fracis_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1:pcnst)



   icwmrdp  => pbuf(icwmrdp_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   rprddp  => pbuf(rprddp_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

   icwmrsh  => pbuf(icwmrsh_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   rprdsh  => pbuf(rprdsh_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

! Dec.29.2009. Sungsu
   sh_frac => pbuf(sh_frac_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   dp_frac => pbuf(dp_frac_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   evapcsh  => pbuf(nevapr_shcu_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   evapcdp  => pbuf(nevapr_dpcu_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

   cldc(:ncol,:)  = dp_frac(:ncol,:) + sh_frac(:ncol,:) ! Sungsu included this.
   evapc(:ncol,:) = evapcsh(:ncol,:) + evapcdp(:ncol,:) ! Sungsu included this.
   clds(:ncol,:)  = cldn(:ncol,:) - cldc(:ncol,:)       ! Stratiform cloud fraction


   ! sum deep and shallow convection contributions
   if (cam_physpkg_is('cam5')) then
      ! Dec.29.2009. Sungsu
      conicw(:ncol,:) = (icwmrdp(:ncol,:)*dp_frac(:ncol,:) + icwmrsh(:ncol,:)*sh_frac(:ncol,:))/ &
                        max(0.01_r8, sh_frac(:ncol,:) + dp_frac(:ncol,:))
   else
      conicw(:ncol,:) = icwmrdp(:ncol,:) + icwmrsh(:ncol,:)
   end if

   cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)


!   fields needed for wet scavenging
    call clddiag( state%t, state%pmid, state%pdel, cmfdqr, evapc, cldn, cldc, clds, cme, evapr, prain, &
         cldv, cldvcu, cldvst, rainmr, ncol )

    ptend%name = 'wetdep'

#if ( defined MODAL_AERO )
    if ( associated(pbuf(dgnum_idx)%fld_ptr) ) then
       dgnum_pbuf => pbuf(dgnum_idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1:ntot_amode )
    else
       call endrun( 'pbuf for DGNUM not allocated in aerosol_wet_intr' )
    end if
    if ( associated(pbuf(dgnumwet_idx)%fld_ptr) ) then
       dgnumwet_pbuf => pbuf(dgnumwet_idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1:ntot_amode )
    else
       call endrun( 'pbuf for DGNUMWET not allocated in aerosol_wet_intr' )
    end if
    if ( associated(pbuf(wetdens_ap_idx)%fld_ptr) ) then
       wetdens_pbuf => pbuf(wetdens_ap_idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1:ntot_amode )
    else
       call endrun( 'pbuf for WETDENS_AP not allocated in aerosol_wet_intr' )
    end if
    if ( associated(pbuf(qaerwat_idx)%fld_ptr) ) then
       qaerwat => pbuf(qaerwat_idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1:ntot_amode )
    else
       call endrun( 'pbuf for QAERWAT not allocated in aerosol_wet_intr' )
    end if
    rate1ord_cw2pr_st => pbuf(rate1_cw2pr_st_idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1 )  ! rce 2010/05/01

    call t_startf('calcsize')

    call modal_aero_calcsize_sub(                                &
                       lchnk,   ncol,    nstep,                  &
                       1,                0,                      &
                       .true.,           dt,                     &
                       state%t,          state%pmid,             &
                       state%pdel,       state%q,                &
                       ptend%q,          ptend%lq,               &
                                         dgnum_pbuf              )

    call t_stopf('calcsize')
    call t_startf('wateruptake')

    call modal_aero_wateruptake_sub(                             &
                       lchnk,   ncol,    nstep,                  &
                       1,                0,                      &
                       .true.,           .true.,                 &
                       dt,               state%q(:,:,1),         &
                       state%t,          state%pmid,             &
                       state%pdel, cldn, state%q,                &
                       ptend%q,          ptend%lq,               &
                       qaerwat,                                  &
                       dgnum_pbuf,       dgnumwet_pbuf,          &
                       wetdens_pbuf                              )

    call t_stopf('wateruptake')
    call mz_aero_wet_intr (state, ptend, nstep, dt, cme, prain, &
            evapr, cldv, cldvcu, cldvst, cldc, cldn, fracis, calday, cmfdqr, evapc, conicw, rainmr, &
            rate1ord_cw2pr_st, &   ! rce 2010/05/01
            dgnumwet_pbuf, qaerwat, cam_out, dlf)


#else
    if ( dust ) then
       !   wet scavenging for dust
       call dust_wet_intr (state, ptend, nstep, dt, lat, clat, cme, prain, &
            evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr, cam_out)
    endif

    if ( progsslt ) then
       !   wet scavenging for prognostic sea salts
       call progseasalts_wet_intr (state, ptend, nstep, dt, lat, clat, cme, prain, &
            evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr)
    endif

!   wet scavenging for mozart aerosols
! can not be done under trop_mozart chem_timestep_tend --
! this need to be done before deep convection so that fracis is set correctly
! for the mozart aerosols -- fracis is used in deep convection routine
    call mz_aero_wet_intr (state, ptend, nstep, dt, cme, prain, &
            evapr, cldv, cldvcu, cldvst, cldc, cldn, fracis, calday, cmfdqr, evapc, conicw, rainmr, cam_out, dlf)
#endif
    return

  end subroutine aerosol_wet_intr

  subroutine aerosol_drydep_intr (state, ptend, cam_in, cam_out, dt, &
                                  fsds, obklen, ustar, prect, pblh, pbuf )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to dry deposition parameterizions and sedimentation of all aerosols
! 
! Method: 
! Use prescribed dry deposition velocities for sulfur and carbon
! Use calculated dry dep velocities from CLM for dust and prognostic seasalt
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use cam_history,       only: outfld
    use physics_types,     only: physics_state, physics_ptend
    use phys_grid,         only: get_lat_all_p, get_rlat_all_p
    use dust_intr,         only: dust_drydep_intr
    use progseasalts_intr, only: progseasalts_drydep_intr

    use time_manager,      only: get_curr_date, get_perp_date, get_curr_calday, &
                                 is_perpetual
    use scyc,              only: scyc_idx1
#if ( defined MODAL_AERO )
    use modal_aero_data,   only: ntot_amode
    use mz_aerosols_intr,  only: mz_aero_dry_intr
#endif
    use phys_buffer,       only: pbuf_size_max, pbuf_fld, pbuf_get_fld_idx

!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)         :: dt      ! time step
    type(physics_state), intent(in )        :: state   ! Physics state variables
    type(cam_in_t),      intent(in), target :: cam_in  ! import state
    type(cam_out_t),     intent(inout)      :: cam_out ! export state
    real(r8), intent(in) :: fsds(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: obklen(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel
    real(r8), intent(in) :: prect(pcols)                 ! prect
    real(r8), intent(in) :: pblh(pcols)                  ! pbl height

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
    type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf

!
! Local
! 
    ! local names of the fields in cam_in
    real(r8), pointer, dimension(:)  :: wvflx	 ! (pcols) water vapor flux
    real(r8), pointer, dimension(:)  :: ts	 ! (pcols) sfc temp
    real(r8), pointer, dimension(:)  :: snowh	 ! (pcols) snow depth
    real(r8), pointer, dimension(:)  :: hflx	 ! (pcols) sensible heat flux
    real(r8), pointer, dimension(:)  :: landfrac ! (pcols) land fraction
    real(r8), pointer, dimension(:)  :: icefrac	 ! (pcols) ice fraction
    real(r8), pointer, dimension(:)  :: ocnfrac	 ! (pcols) ocn fraction


    integer  lat(pcols)                  ! latitude index for S->N storage
    real(r8) clat(pcols)                 ! latitude 
    integer lchnk
    integer ncol
    integer ix
    integer m

    integer :: yr, mon, day, ncsec
    integer :: ncdate

#if ( defined MODAL_AERO )
    real(r8), pointer, dimension(:,:,:) :: dgnumwet_pbuf, wetdens_pbuf
    real(r8), pointer, dimension(:,:,:) :: qaerwat  ! aerosol water
#endif


    lchnk = state%lchnk
    ncol = state%ncol
    
    wvflx    => cam_in%cflx(:,1)
    ts       => cam_in%tref(:)
    snowh    => cam_in%snowhland(:)
    hflx     => cam_in%shf(:)
    landfrac => cam_in%landfrac(:)
    icefrac  => cam_in%icefrac(:)
    ocnfrac  => cam_in%ocnfrac(:)

    if ( is_perpetual() ) then
       call get_perp_date(yr, mon, day, ncsec)
    else
       call get_curr_date(yr, mon, day, ncsec)
    end if

    call get_lat_all_p(lchnk, ncol, lat)
    call get_rlat_all_p(lchnk, ncol, clat)

!   note that tendencies are not only in sfc layer (because of sedimentation)
!   and that ptend is updated within each subroutine for different species

    ptend%name = 'drydep'

#if ( defined MODAL_AERO )
    if ( associated(pbuf(dgnumwet_idx)%fld_ptr) ) then
       dgnumwet_pbuf => pbuf(dgnumwet_idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1:ntot_amode )
    else
       call endrun( 'pbuf for DGNUMWET not allocated in aerosol_drydep_intr' )
    end if
    if ( associated(pbuf(wetdens_ap_idx)%fld_ptr) ) then
       wetdens_pbuf => pbuf(wetdens_ap_idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1:ntot_amode )
    else
       call endrun( 'pbuf for WETDENS_AP not allocated in aerosol_drydep_intr' )
    end if
    if ( associated(pbuf(qaerwat_idx)%fld_ptr) ) then
       qaerwat => pbuf(qaerwat_idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1:ntot_amode )
    else
       call endrun( 'pbuf for QAERWAT not allocated in aerosol_drydep_intr' )
    end if

    call mz_aero_dry_intr (state, ptend, wvflx, dt, lat, clat, &
            fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, mon, &
            landfrac, icefrac, ocnfrac, cam_in%fv, cam_in%ram1, &
            dgnumwet_pbuf, wetdens_pbuf, qaerwat, cam_out)
#else
    if ( dust ) then
       call dust_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
            fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, mon, &
            landfrac, icefrac, ocnfrac, cam_in%fv,  cam_in%ram1, cam_out)
    endif

    if ( progsslt ) then
       call progseasalts_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
            fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, mon, landfrac, &
            icefrac, ocnfrac, cam_in%fv,  cam_in%ram1)
    endif
#endif


    return
  end subroutine aerosol_drydep_intr

#if ( defined MODAL_AERO )
  subroutine aerosol_emis_intr (state, ptend, cflx, dt,ocnfrc,sst)
#else
  subroutine aerosol_emis_intr (state, ptend, cflx, dt,ocnfrc)
#endif

!----------------------------------------------------------------------- 
! 
! Purpose: 
! return surface fluxes of aerosol species and tendencies in surface layer
! due to surface emissions
! 
! Method: 
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use cam_history,        only: outfld
    use physics_types,      only: physics_state, physics_ptend
    use phys_grid,          only: get_lon_all_p, get_lat_all_p, get_rlat_all_p
    use time_manager,       only: get_curr_date, get_perp_date, get_curr_calday, &
                                  is_perpetual
    use dust_intr,          only: dust_emis_intr
    use progseasalts_intr,  only: progseasalts_emis_intr

!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables

     real(r8),            intent(in)  :: ocnfrc(pcols)

#if ( defined MODAL_AERO )
    real(r8),            intent(in)  :: sst(pcols)     ! Sea surface temperature
#endif
    type(physics_ptend), intent(inout) :: ptend        ! indivdual parameterization tendencies
    real(r8), intent(inout)  :: cflx(pcols,pcnst)      ! Surface constituent flux (kg/m^2/s)

    integer  lat(pcols)                  ! latitude index 
    integer  lon(pcols)                  ! longitude index
    real(r8) clat(pcols)                 ! latitude 
    integer lchnk
    integer ncol
    integer i
    integer m
    real(r8) astmp(pcols,pver,3)
!
    real(r8) :: calday        ! current calendar day
    integer :: yr, mon, day, ncsec
    integer :: ncdate
    real(r8) :: so2sf(pcols), so4sf(pcols), dmssf(pcols)
    integer :: k

    ptend%name = 'aerosol_emis'

#if ( defined MODAL_AERO )
    call dust_emis_intr(state,ptend,dt,cflx)
    call progseasalts_emis_intr(state,ptend,dt,cflx,ocnfrc,sst)
#else
    if ( dust ) then
       call dust_emis_intr(state,ptend,dt,cflx)
    endif
    if ( progsslt ) then
       call progseasalts_emis_intr(state,ptend,dt,cflx,ocnfrc)
    endif
#endif


    return
  end subroutine aerosol_emis_intr 

end module aerosol_intr
