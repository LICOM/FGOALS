  module macrop_driver

  !-------------------------------------------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the CAM interface to the prognostic cloud macrophysics
  !
  ! Author: Andrew Gettelman, Cheryl Craig October 2010
  ! Origin: modified from stratiform.F90 elements 
  !    (Boville 2002, Coleman 2004, Park 2009, Kay 2010)
  !-------------------------------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: latice
  use phys_control,  only: phys_getopts
  use constituents,  only: cnst_get_ind
  use perf_mod,      only: t_startf, t_stopf
  use cam_logfile,   only: iulog

  implicit none
  private
  save

  public :: macrop_driver_register
  public :: macrop_driver_init
  public :: macrop_driver_tend

  ! ------------------------- !
  ! Private Module Parameters !
  ! ------------------------- !

  ! 'cu_det_st' : If .true. (.false.), detrain cumulus liquid condensate into the pre-existing liquid stratus 
  !               (environment) without (with) macrophysical evaporation. If there is no pre-esisting stratus, 
  !               evaporate cumulus liquid condensate. This option only influences the treatment of cumulus
  !               liquid condensate, not cumulus ice condensate.

    logical,          private, parameter :: cu_det_st  = .false.  

  ! -------------------------------- !
  ! End of Private Module Parameters !
  ! -------------------------------- !

  logical            :: use_shfrc                       ! Local copy of flag from convect_shallow_use_shfrc
  character(len=16)  :: microp_scheme                   ! Microphysics scheme

  integer :: &
    ixcldliq,     &! cloud liquid amount index
    ixcldice,     &! cloud ice amount index
    ixnumliq,     &! cloud liquid number index
    ixnumice,     &! cloud ice water index
    qcwat_idx,    &! qcwat index in physics buffer
    lcwat_idx,    &! lcwat index in physics buffer
    iccwat_idx,   &! iccwat index in physics buffer
    nlwat_idx,    &! nlwat index in physics buffer
    niwat_idx,    &! niwat index in physics buffer
    tcwat_idx,    &! tcwat index in physics buffer
    CC_T_idx,     &!
    CC_qv_idx,    &!
    CC_ql_idx,    &!
    CC_qi_idx,    &!
    CC_nl_idx,    &!
    CC_ni_idx,    &!
    CC_qlst_idx,  &!
    cld_idx,      &! cld index in physics buffer
    ast_idx,      &! stratiform cloud fraction index in physics buffer
    aist_idx,     &! ice stratiform cloud fraction index in physics buffer
    alst_idx,     &! liquid stratiform cloud fraction index in physics buffer
    qist_idx,     &! ice stratiform in-cloud IWC 
    qlst_idx,     &! liquid stratiform in-cloud LWC  
    concld_idx,   &! concld index in physics buffer
    rhdfda_idx,   &! rhdfda index in physics buffer
    fice_idx,     &  
    qme_idx,      & 
    shfrc_idx,    & 
    rhu00_idx,    & 
    concldql_idx, & 
    sh_frac_idx,  &
    dp_frac_idx,  &
    qini_idx,     &
    cldliqini_idx,&
    cldiceini_idx,&
    tini_idx      
    

  contains

  ! ===============================================================================

  subroutine macrop_driver_register

  !---------------------------------------------------------------------- !
  !                                                                       !
  ! Register the constituents (cloud liquid and cloud ice) and the fields !
  ! in the physics buffer.                                                !
  !                                                                       !
  !---------------------------------------------------------------------- !

   use phys_buffer,  only: pbuf_times, pbuf_add

  !-----------------------------------------------------------------------
    call phys_getopts(microp_scheme_out=microp_scheme)

    call pbuf_add('AST',     'global',  1, pver, pbuf_times,     ast_idx)
    call pbuf_add('AIST',    'global',  1, pver, pbuf_times,    aist_idx)
    call pbuf_add('ALST',    'global',  1, pver, pbuf_times,    alst_idx)
    call pbuf_add('QIST',    'global',  1, pver, pbuf_times,    qist_idx)
    call pbuf_add('QLST',    'global',  1, pver, pbuf_times,    qlst_idx)
    call pbuf_add('CONCLD',  'global',  1, pver, pbuf_times,  concld_idx)
    call pbuf_add('RHDFDA',  'global',  1, pver, pbuf_times,  rhdfda_idx)
    call pbuf_add('RHU00',   'global',  1, pver, pbuf_times,   rhu00_idx)
    call pbuf_add('QCWAT',   'global',  1, pver, pbuf_times,   qcwat_idx)
    call pbuf_add('LCWAT',   'global',  1, pver, pbuf_times,   lcwat_idx)
    call pbuf_add('ICCWAT',  'global',  1, pver, pbuf_times,  iccwat_idx)
    call pbuf_add('NLWAT',   'global',  1, pver, pbuf_times,   nlwat_idx)
    call pbuf_add('NIWAT',   'global',  1, pver, pbuf_times,   niwat_idx)
    call pbuf_add('TCWAT',   'global',  1, pver, pbuf_times,   tcwat_idx)
    call pbuf_add('CC_T',    'global',  1, pver, pbuf_times,    CC_T_idx)
    call pbuf_add('CC_qv',   'global',  1, pver, pbuf_times,   CC_qv_idx)
    call pbuf_add('CC_ql',   'global',  1, pver, pbuf_times,   CC_ql_idx)
    call pbuf_add('CC_qi',   'global',  1, pver, pbuf_times,   CC_qi_idx)
    call pbuf_add('CC_nl',   'global',  1, pver, pbuf_times,   CC_nl_idx)
    call pbuf_add('CC_ni',   'global',  1, pver, pbuf_times,   CC_ni_idx)
    call pbuf_add('CC_qlst', 'global',  1, pver, pbuf_times, CC_qlst_idx)
    call pbuf_add('CLD',     'global',  1, pver, pbuf_times,     cld_idx)
    call pbuf_add('CONCLDQL','physpkg', 1, pver, 1,         concldql_idx) 

    call pbuf_add('SH_FRAC', 'physpkg', 1, pver, 1,          sh_frac_idx)
    call pbuf_add('DP_FRAC', 'physpkg', 1, pver, 1,          dp_frac_idx)

    call pbuf_add('QINI',      'physpkg',1,pver, 1,             qini_idx)
    call pbuf_add('CLDLIQINI', 'physpkg',1,pver, 1,        cldliqini_idx)
    call pbuf_add('CLDICEINI', 'physpkg',1,pver, 1,        cldiceini_idx)
    call pbuf_add('TINI',      'physpkg',1,pver, 1,             tini_idx)


  end subroutine macrop_driver_register

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine macrop_driver_init

  !-------------------------------------------- !
  !                                             !
  ! Initialize the cloud water parameterization !
  !                                             ! 
  !-------------------------------------------- !

    use cam_history,     only: addfld, add_default, phys_decomp
    use convect_shallow, only: convect_shallow_use_shfrc
    use phys_buffer,     only: pbuf_get_fld_idx


    logical              :: history_aerosol      ! Output the MAM aerosol tendencies
    logical              :: history_microphysics ! Output variables for microphysics diagnostics package
    logical              :: history_budget       ! Output tendencies and state variables for CAM4
                                                 ! temperature, water vapor, cloud ice and cloud
                                                 ! liquid budgets.
    integer              :: history_budget_histfile_num ! output history file number for budget fields

  !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out        = history_aerosol      , &
                       history_microphysics_out   = history_microphysics , & 
                       history_budget_out         = history_budget       , &
                       history_budget_histfile_num_out = history_budget_histfile_num)

  ! Initialization routine for cloud macrophysics

  ! Find out whether shfrc from convect_shallow will be used in cldfrc

    if( convect_shallow_use_shfrc() ) then
        use_shfrc = .true.
        shfrc_idx   = pbuf_get_fld_idx('shfrc')
   else 
        use_shfrc = .false.
    endif


    call addfld ('DPDLFLIQ ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from deep convection'             ,phys_decomp)
    call addfld ('DPDLFICE ', 'kg/kg/s ', pver, 'A', 'Detrained ice from deep convection'                      ,phys_decomp)
    call addfld ('SHDLFLIQ ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from shallow convection'          ,phys_decomp)
    call addfld ('SHDLFICE ', 'kg/kg/s ', pver, 'A', 'Detrained ice from shallow convection'                   ,phys_decomp)
    call addfld ('DPDLFT   ', 'K/s     ', pver, 'A', 'T-tendency due to deep convective detrainment'           ,phys_decomp)
    call addfld ('SHDLFT   ', 'K/s     ', pver, 'A', 'T-tendency due to shallow convective detrainment'        ,phys_decomp)

    call addfld ('ZMDLF    ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from ZM convection'               ,phys_decomp)

    call addfld ('MACPDT   ', 'W/kg    ', pver, 'A', 'Heating tendency - Revised  macrophysics'                ,phys_decomp)
    call addfld ('MACPDQ   ', 'kg/kg/s ', pver, 'A', 'Q tendency - Revised macrophysics'                       ,phys_decomp)
    call addfld ('MACPDLIQ ', 'kg/kg/s ', pver, 'A', 'CLDLIQ tendency - Revised macrophysics'                  ,phys_decomp)
    call addfld ('MACPDICE ', 'kg/kg/s ', pver, 'A', 'CLDICE tendency - Revised macrophysics'                  ,phys_decomp)

    call addfld ('CLDVAPADJ', 'kg/kg/s ', pver, 'A', 'Q tendency associated with liq/ice adjustment - Revised macrophysics' ,phys_decomp)
    call addfld ('CLDLIQADJ', 'kg/kg/s ', pver, 'A', 'CLDLIQ adjustment tendency - Revised macrophysics'       ,phys_decomp)
    call addfld ('CLDICEADJ', 'kg/kg/s ', pver, 'A', 'CLDICE adjustment tendency - Revised macrophysics'       ,phys_decomp)
    call addfld ('CLDLIQDET', 'kg/kg/s ', pver, 'A', 'Detrainment of conv cld liq into envrionment  - Revised macrophysics' ,phys_decomp)
    call addfld ('CLDICEDET', 'kg/kg/s ', pver, 'A', 'Detrainment of conv cld ice into envrionment  - Revised macrophysics' ,phys_decomp)
    call addfld ('CLDLIQLIM', 'kg/kg/s ', pver, 'A', 'CLDLIQ limiting tendency - Revised macrophysics'         ,phys_decomp)
    call addfld ('CLDICELIM', 'kg/kg/s ', pver, 'A', 'CLDICE limiting tendency - Revised macrophysics'         ,phys_decomp)

    call addfld ('CNVCLD   ', 'fraction', 1,    'A', 'Vertically integrated convective cloud amount'           ,phys_decomp)
    call addfld ('CLDST    ', 'fraction', pver, 'A', 'Stratus cloud fraction'                                  ,phys_decomp)
    call addfld ('CONCLD   ', 'fraction', pver, 'A', 'Convective cloud cover'                                  ,phys_decomp)
 
    call addfld ('CLDLIQSTR   ', 'kg/kg', pver, 'A', 'Stratiform CLDLIQ'                                  ,phys_decomp)
    call addfld ('CLDICESTR   ', 'kg/kg', pver, 'A', 'Stratiform CLDICE'                                  ,phys_decomp)
    call addfld ('CLDLIQCON   ', 'kg/kg', pver, 'A', 'Convective CLDLIQ'                                  ,phys_decomp)
    call addfld ('CLDICECON   ', 'kg/kg', pver, 'A', 'Convective CLDICE'                                  ,phys_decomp)

    if ( history_budget ) then

          call add_default ('DPDLFLIQ ', history_budget_histfile_num, ' ')
          call add_default ('DPDLFICE ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFLIQ ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFICE ', history_budget_histfile_num, ' ')
          call add_default ('DPDLFT   ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFT   ', history_budget_histfile_num, ' ')
          call add_default ('ZMDLF    ', history_budget_histfile_num, ' ')

          call add_default ('MACPDT   ', history_budget_histfile_num, ' ')
          call add_default ('MACPDQ   ', history_budget_histfile_num, ' ')
          call add_default ('MACPDLIQ ', history_budget_histfile_num, ' ')
          call add_default ('MACPDICE ', history_budget_histfile_num, ' ')
 
          call add_default ('CLDVAPADJ', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQLIM', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQDET', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQADJ', history_budget_histfile_num, ' ')
          call add_default ('CLDICELIM', history_budget_histfile_num, ' ')
          call add_default ('CLDICEDET', history_budget_histfile_num, ' ')
          call add_default ('CLDICEADJ', history_budget_histfile_num, ' ')


    end if

! Retrieve the indices for pbuf entries

    fice_idx    = pbuf_get_fld_idx('FICE')
    qme_idx     = pbuf_get_fld_idx('QME')

    return
  end subroutine macrop_driver_init

  !============================================================================ !
  !                                                                             !
  !============================================================================ !


  subroutine macrop_driver_tend(                             &
             state, ptend_all, dtime, landfrac,  &
             ocnfrac,  snowh,                       &
             dlf, dlf2, cmfmc, cmfmc2, ts,          &
             sst, zdu,       &
             pbuf, state_eq,cmeliq)

  !-------------------------------------------------------- !  
  !                                                         ! 
  ! Purpose:                                                !
  !                                                         !
  ! Interface to detrain, cloud fraction and                !
  !     cloud macrophysics subroutines                      !
  !                                                         ! 
  ! Author: A. Gettelman, C. Craig, Oct 2010                !
  ! based on stratiform_tend by D.B. Coleman 4/2010         !
  !                                                         !
  !-------------------------------------------------------- !

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use cloud_fraction,   only: cldfrc
  use physics_types,    only: physics_state, physics_ptend, physics_tend
  use physics_types,    only: physics_ptend_init, physics_update, physics_tend_init
  use physics_types,    only: physics_ptend_sum,  physics_state_copy
  use cam_history,      only: outfld
  use phys_buffer,      only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx
  use constituents,     only: cnst_get_ind, pcnst
  use cldwat,           only: pcond, cldwat_fice
  use cldwat2m_macro,   only: mmacro_pcond
  use physconst,        only: cpair
  use time_manager,     only: get_nstep

  implicit none

  !
  ! Input arguments
  !

  type(physics_state), intent(in)    :: state       ! State variables
  type(physics_ptend), intent(out)   :: ptend_all   ! Package tendencies
  type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf

  real(r8), intent(in)  :: dtime                    ! Timestep
  real(r8), intent(in)  :: landfrac(pcols)          ! Land fraction (fraction)
  real(r8), intent(in)  :: ocnfrac (pcols)          ! Ocean fraction (fraction)
  real(r8), intent(in)  :: snowh(pcols)             ! Snow depth over land, water equivalent (m)
  real(r8), intent(in)  :: dlf(pcols,pver)          ! Detrained water from convection schemes
  real(r8), intent(in)  :: dlf2(pcols,pver)         ! Detrained water from shallow convection scheme
  real(r8), intent(in)  :: cmfmc(pcols,pverp)       ! Deep + Shallow Convective mass flux [ kg /s/m^2 ]
  real(r8), intent(in)  :: cmfmc2(pcols,pverp)      ! Shallow convective mass flux [ kg/s/m^2 ]

  real(r8), intent(in)  :: ts(pcols)                ! Surface temperature
  real(r8), intent(in)  :: sst(pcols)               ! Sea surface temperature
  real(r8), intent(in)  :: zdu(pcols,pver)          ! Detrainment rate from deep convection

  ! Equilibrium state variables at the end of macrophysics
  ! Below 'state_eq' is for future use as the input of radiation'PBL scheme

  type(physics_state), intent(out) :: state_eq   ! Equilibrium state variables at the end of macrophysics

  ! for passing to microphysics
  real(r8), intent(out) :: cmeliq(pcols,pver)                      ! Rate of cond-evap of liq within the cloud


  !
  ! Local variables
  !

  type(physics_state)   :: state1                   ! Local copy of the state variable
  type(physics_tend )   :: tend                     ! Physics tendencies (empty, needed for physics_update call)
  type(physics_ptend)   :: ptend_loc                ! Package tendencies

  integer i,k
  integer :: lchnk                                  ! Chunk identifier
  integer :: ncol                                   ! Number of atmospheric columns
  integer :: conv_water_in_rad

  ! Physics buffer fields

  integer itim
  real(r8), pointer, dimension(:,:) :: rhdfda       !
  real(r8), pointer, dimension(:,:) :: rhu00        ! 
  real(r8), pointer, dimension(:,:) :: qcwat        ! Cloud water old q
  real(r8), pointer, dimension(:,:) :: tcwat        ! Cloud water old temperature
  real(r8), pointer, dimension(:,:) :: lcwat        ! Cloud liquid water old q
  real(r8), pointer, dimension(:,:) :: iccwat       ! Cloud ice water old q
  real(r8), pointer, dimension(:,:) :: nlwat        ! Cloud liquid droplet number condentration. old.
  real(r8), pointer, dimension(:,:) :: niwat        ! Cloud ice    droplet number condentration. old.
  real(r8), pointer, dimension(:,:) :: CC_T         ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qv        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_ql        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qi        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_nl        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_ni        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qlst      ! In-liquid stratus microphysical tendency
  real(r8), pointer, dimension(:,:) :: cld          ! Total cloud fraction
  real(r8), pointer, dimension(:,:) :: ast          ! Relative humidity cloud fraction
  real(r8), pointer, dimension(:,:) :: aist         ! Physical ice stratus fraction
  real(r8), pointer, dimension(:,:) :: alst         ! Physical liquid stratus fraction
  real(r8), pointer, dimension(:,:) :: qist         ! Physical in-cloud IWC
  real(r8), pointer, dimension(:,:) :: qlst         ! Physical in-cloud LWC
  real(r8), pointer, dimension(:,:) :: concld       ! Convective cloud fraction
  real(r8), pointer, dimension(:,:) :: qme

  real(r8) :: shfrc(pcols,pver)                     ! Cloud fraction from shallow convection scheme

  ! Convective cloud to the physics buffer for purposes of ql contrib. to radn.

  real(r8), pointer, dimension(:,:) :: concld_ql    ! Convective cloud
  real(r8), pointer, dimension(:,:) :: fice_ql      ! Cloud ice/water partitioning ratio.

  ! Local variables for cldfrc

  real(r8)  cldst(pcols,pver)                       ! Stratus cloud fraction
  real(r8)  rhcloud(pcols,pver)                     ! Relative humidity cloud (last timestep)
  real(r8)  rhcloud2(pcols,pver)                    ! Relative humidity cloud (perturbation)
  real(r8)  clc(pcols)                              ! Column convective cloud amount
  real(r8)  relhum(pcols,pver)                      ! RH, output to determine drh/da
  real(r8)  rhu002(pcols,pver)                      ! Same as rhu00 but for perturbed rh 
  real(r8)  cld2(pcols,pver)                        ! Same as cld but for perturbed rh
  real(r8)  concld2(pcols,pver)                     ! Same as concld but for perturbed rh 
  real(r8)  cldst2(pcols,pver)                      ! Same as cldst but for perturbed rh 
  real(r8)  relhum2(pcols,pver)                     ! RH after  perturbation            
  real(r8)  icecldf(pcols,pver)                     ! Ice cloud fraction
  real(r8)  liqcldf(pcols,pver)                     ! Liquid cloud fraction (combined into cloud)
  real(r8)  icecldf2(pcols,pver)                    ! Ice cloud fraction
  real(r8)  liqcldf2(pcols,pver)                    ! Liquid cloud fraction (combined into cloud)

  ! Local variables for macrophysics

  real(r8)  rdtime                                  ! 1./dtime
  real(r8)  qtend(pcols,pver)                       ! Moisture tendencies
  real(r8)  ttend(pcols,pver)                       ! Temperature tendencies
  real(r8)  ltend(pcols,pver)                       ! Cloud liquid water tendencies
  real(r8)  fice(pcols,pver)                        ! Fractional ice content within cloud
  real(r8)  fsnow(pcols,pver)                       ! Fractional snow production
  real(r8)  homoo(pcols,pver)  
  real(r8)  qcreso(pcols,pver)  
  real(r8)  prcio(pcols,pver)  
  real(r8)  praio(pcols,pver)  
  real(r8)  qireso(pcols,pver)
  real(r8)  ftem(pcols,pver)
  real(r8)  pracso (pcols,pver) 
  real(r8)  dpdlfliq(pcols,pver)
  real(r8)  dpdlfice(pcols,pver)
  real(r8)  shdlfliq(pcols,pver)
  real(r8)  shdlfice(pcols,pver)
  real(r8)  dpdlft  (pcols,pver)
  real(r8)  shdlft  (pcols,pver)

  real(r8)  dum1
  real(r8)  qc(pcols,pver)
  real(r8)  qi(pcols,pver)
  real(r8)  nc(pcols,pver)
  real(r8)  ni(pcols,pver)

  ! Output from mmacro_pcond

  real(r8)  tlat(pcols,pver)
  real(r8)  qvlat(pcols,pver)
  real(r8)  qcten(pcols,pver)
  real(r8)  qiten(pcols,pver)
  real(r8)  ncten(pcols,pver)
  real(r8)  niten(pcols,pver)

  ! Output from mmacro_pcond

  real(r8)  qvadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (vapor)
  real(r8)  qladj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (liquid)
  real(r8)  qiadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (ice)
  real(r8)  qllim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (liquid)
  real(r8)  qilim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (ice)

  ! For revised macophysics, mmacro_pcond

  real(r8)  itend(pcols,pver)
  real(r8)  lmitend(pcols,pver)
  real(r8)  zeros(pcols,pver)
  real(r8)  t_inout(pcols,pver)
  real(r8)  qv_inout(pcols,pver)
  real(r8)  ql_inout(pcols,pver)
  real(r8)  qi_inout(pcols,pver)
  real(r8)  concld_old(pcols,pver)

  real(r8)  nl_inout(pcols,pver)
  real(r8)  ni_inout(pcols,pver)

  real(r8)  nltend(pcols,pver)
  real(r8)  nitend(pcols,pver)


  ! For detraining cumulus condensate into the 'stratus' without evaporation
  ! This is for use in mmacro_pcond

  real(r8)  dlf_T(pcols,pver)
  real(r8)  dlf_qv(pcols,pver)
  real(r8)  dlf_ql(pcols,pver)
  real(r8)  dlf_qi(pcols,pver)
  real(r8)  dlf_nl(pcols,pver)
  real(r8)  dlf_ni(pcols,pver)

  ! Local variables for CFMIP calculations
  real(r8) :: mr_lsliq(pcols,pver)			! mixing_ratio_large_scale_cloud_liquid (kg/kg)
  real(r8) :: mr_lsice(pcols,pver)			! mixing_ratio_large_scale_cloud_ice (kg/kg)
  real(r8) :: mr_ccliq(pcols,pver)			! mixing_ratio_convective_cloud_liquid (kg/kg)
  real(r8) :: mr_ccice(pcols,pver)			! mixing_ratio_convective_cloud_ice (kg/kg)

  ! ======================================================================

  lchnk = state%lchnk
  ncol  = state%ncol

  call cnst_get_ind('CLDLIQ', ixcldliq)
  call cnst_get_ind('CLDICE', ixcldice)
  call cnst_get_ind('NUMLIQ', ixnumliq)
  call cnst_get_ind('NUMICE', ixnumice)

  call phys_getopts( conv_water_in_rad_out = conv_water_in_rad )

  call physics_state_copy(state,state1)             ! Copy state to local state1.
  call physics_ptend_init(ptend_loc)                ! Initialize local ptend type
  call physics_ptend_init(ptend_all)                ! Initialize output ptend type
  call physics_tend_init(tend)                      ! tend here is just a null place holder

  ! Associate pointers with physics buffer fields

  itim = pbuf_old_tim_idx()

  qcwat => pbuf(qcwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  rhdfda => pbuf(rhdfda_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  rhu00 => pbuf(rhu00_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  tcwat => pbuf(tcwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  lcwat => pbuf(lcwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  iccwat => pbuf(iccwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  nlwat => pbuf(nlwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  niwat => pbuf(niwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_T => pbuf(cc_t_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_qv => pbuf(cc_qv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_ql => pbuf(cc_ql_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_qi => pbuf(cc_qi_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_nl => pbuf(cc_nl_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_ni => pbuf(cc_ni_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_qlst => pbuf(cc_qlst_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  cld => pbuf(cld_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ast => pbuf(ast_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  aist => pbuf(aist_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  alst => pbuf(alst_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  qist => pbuf(qist_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  qlst => pbuf(qlst_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  concld => pbuf(concld_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

! For purposes of convective ql.

  concld_ql => pbuf(concldql_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)                                      

  fice_ql => pbuf(fice_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

  qme  => pbuf(qme_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)


  ! Initialize convective detrainment tendency

  dlf_T(:,:)  = 0._r8
  dlf_qv(:,:) = 0._r8
  dlf_ql(:,:) = 0._r8
  dlf_qi(:,:) = 0._r8
  dlf_nl(:,:) = 0._r8
  dlf_ni(:,:) = 0._r8

   ! ------------------------------------- !
   ! From here, process computation begins ! 
   ! ------------------------------------- !

   ! ----------------------------------------------------------------------------- !
   ! Detrainment of convective condensate into the environment or stratiform cloud !
   ! ----------------------------------------------------------------------------- !

   ptend_loc%name = 'pcwdetrain'

   ptend_loc%lq(ixcldliq) = .TRUE.
   ptend_loc%lq(ixcldice) = .TRUE.
   ptend_loc%lq(ixnumliq) = .TRUE.
   ptend_loc%lq(ixnumice) = .TRUE.
   ptend_loc%ls           = .TRUE.

     ! Procedures :
     ! (1) Partition detrained convective cloud water into liquid and ice based on T.
     !     This also involves heating.
     !     If convection scheme can handle this internally, this step is not necssary.
     ! (2) Assuming a certain effective droplet radius, computes number concentration
     !     of detrained convective cloud liquid and ice.
     ! (3) If 'cu_det_st = .true' ('false'), detrain convective cloud 'liquid' into 
     !     the pre-existing 'liquid' stratus ( mean environment ).  The former does
     !     not involve any macrophysical evaporation while the latter does. This is
     !     a kind of 'targetted' deposition. Then, force in-stratus LWC to be bounded 
     !     by qcst_min and qcst_max in mmacro_pcond.
     ! (4) In contrast to liquid, convective ice is detrained into the environment 
     !     and involved in the sublimation. Similar bounds as liquid stratus are imposed.
     ! This is the key procesure generating upper-level cirrus clouds.
     ! The unit of dlf : [ kg/kg/s ]

   do k = 1, pver
   do i = 1, state1%ncol
      if( state1%t(i,k) > 268.15_r8 ) then
          dum1 = 0.0_r8
      elseif( state1%t(i,k) < 238.15_r8 ) then
          dum1 = 1.0_r8
      else
          dum1 = ( 268.15_r8 - state1%t(i,k) ) / 30._r8
      endif
      ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * ( 1._r8 - dum1 )
      ptend_loc%q(i,k,ixcldice) = dlf(i,k) * dum1
    ! dum2                      = dlf(i,k) * ( 1._r8 - dum1 )
      ptend_loc%q(i,k,ixnumliq) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) * ( 1._r8 - dum1 ) ) / (4._r8*3.14_r8* 8.e-6_r8**3*997._r8) + & ! Deep    Convection
                                  3._r8 * (                         dlf2(i,k)    * ( 1._r8 - dum1 ) ) / (4._r8*3.14_r8*10.e-6_r8**3*997._r8)     ! Shallow Convection 
    ! dum2                      = dlf(i,k) * dum1
      ptend_loc%q(i,k,ixnumice) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) *  dum1 ) / (4._r8*3.14_r8*25.e-6_r8**3*500._r8) + & ! Deep    Convection
                                  3._r8 * (                         dlf2(i,k)    *  dum1 ) / (4._r8*3.14_r8*50.e-6_r8**3*500._r8)     ! Shallow Convection
      ptend_loc%s(i,k)          = dlf(i,k) * dum1 * latice
    ! Targetted detrainment of convective liquid water either directly into the
    ! existing liquid stratus or into the environment. 
      if( cu_det_st ) then
          dlf_T(i,k)  = ptend_loc%s(i,k)/cpair
          dlf_qv(i,k) = 0._r8
          dlf_ql(i,k) = ptend_loc%q(i,k,ixcldliq)
          dlf_qi(i,k) = ptend_loc%q(i,k,ixcldice)
          dlf_nl(i,k) = ptend_loc%q(i,k,ixnumliq)
          dlf_ni(i,k) = ptend_loc%q(i,k,ixnumice)         
          ptend_loc%q(i,k,ixcldliq) = 0._r8
          ptend_loc%q(i,k,ixcldice) = 0._r8
          ptend_loc%q(i,k,ixnumliq) = 0._r8
          ptend_loc%q(i,k,ixnumice) = 0._r8
          ptend_loc%s(i,k)          = 0._r8
          dpdlfliq(i,k)             = 0._r8
          dpdlfice(i,k)             = 0._r8
          shdlfliq(i,k)             = 0._r8
          shdlfice(i,k)             = 0._r8
          dpdlft  (i,k)             = 0._r8
          shdlft  (i,k)             = 0._r8
       else
          dpdlfliq(i,k) = ( dlf(i,k) - dlf2(i,k) ) * ( 1._r8 - dum1 )
          dpdlfice(i,k) = ( dlf(i,k) - dlf2(i,k) ) * ( dum1 )
          shdlfliq(i,k) = dlf2(i,k) * ( 1._r8 - dum1 )
          shdlfice(i,k) = dlf2(i,k) * ( dum1 )
          dpdlft  (i,k) = ( dlf(i,k) - dlf2(i,k) ) * dum1 * latice/cpair
          shdlft  (i,k) = dlf2(i,k) * dum1 * latice/cpair
      endif
   end do
   end do

   call outfld( 'DPDLFLIQ ', dpdlfliq, pcols, lchnk )
   call outfld( 'DPDLFICE ', dpdlfice, pcols, lchnk )
   call outfld( 'SHDLFLIQ ', shdlfliq, pcols, lchnk )
   call outfld( 'SHDLFICE ', shdlfice, pcols, lchnk )
   call outfld( 'DPDLFT   ', dpdlft  , pcols, lchnk )
   call outfld( 'SHDLFT   ', shdlft  , pcols, lchnk )

   call outfld( 'ZMDLF',     dlf     , pcols, state1%lchnk )

 ! Add hie detrainment tendency to tend from the other prior processes

   call physics_ptend_sum( ptend_loc, ptend_all, state )
   call physics_update( state1, tend, ptend_loc, dtime )
   call physics_ptend_init( ptend_loc )

   ! -------------------------------------- !
   ! Computation of Various Cloud Fractions !
   ! -------------------------------------- !

   ! ----------------------------------------------------------------------------- !
   ! Treatment of cloud fraction in CAM4 and CAM5 differs                          !  
   ! (1) CAM4                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( empirical fcn of mass flux )        !
   !     . Stratus AMT = max( RH stratus AMT, Stability Stratus AMT )              !
   !     . Cumulus and Stratus are 'minimally' overlapped without hierarchy.       !
   !     . Cumulus LWC,IWC is assumed to be the same as Stratus LWC,IWC            !
   ! (2) CAM5                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( internally fcn of mass flux and w ) !
   !     . Stratus AMT = fcn of environmental-mean RH ( no Stability Stratus )     !
   !     . Cumulus and Stratus are non-overlapped with higher priority on Cumulus  !
   !     . Cumulus ( both Deep and Shallow ) has its own LWC and IWC.              !
   ! ----------------------------------------------------------------------------- ! 

   concld_old(:ncol,:pver) = concld(:ncol,:pver)

   if( use_shfrc ) then
       shfrc(:pcols,:pver) = pbuf(shfrc_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   else 
       shfrc(:,:) = 0._r8
   endif

   ! CAM5 only uses 'concld' output from the below subroutine. 
   ! Stratus ('ast' = max(alst,aist)) and total cloud fraction ('cld = ast + concld')
   ! will be computed using this updated 'concld' in the stratiform macrophysics 
   ! scheme (mmacro_pcond) later below. 
   ! Note 'shfrc' and' deep convective cloud fraction' will be saved into the 
   ! physical buffer (SH_FRAC,DP_FRAC) within cldfrc.

   call t_startf("cldfrc")
   call cldfrc( lchnk, ncol, pbuf,                                                 &
                state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                shfrc, use_shfrc,                                                  &
                cld, rhcloud, clc, state1%pdel,                                    &
                cmfmc, cmfmc2, landfrac,snowh, concld, cldst,                      &
                ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu00,                &
                state1%q(:,:,ixcldice), icecldf, liqcldf,                          &
                relhum, 0 )    

   ! Re-calculate cloud with perturbed rh add call cldfrc  

   call cldfrc( lchnk, ncol, pbuf,                                                 &
                state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                shfrc, use_shfrc,                                                  &
                cld2, rhcloud2, clc, state1%pdel,                                  &
                cmfmc, cmfmc2, landfrac, snowh, concld2, cldst2,                   &
                ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu002,               &
                state1%q(:,:,ixcldice), icecldf2, liqcldf2,                        &
                relhum2, 1 )              

   call t_stopf("cldfrc")

   ! Add following to estimate rhdfda. Below block is only for CAM4                       

   rhu00(:ncol,1) = 2.0_r8

   do k = 1, pver
   do i = 1, ncol
      if( relhum(i,k) < rhu00(i,k) ) then
          rhdfda(i,k) = 0.0_r8
      elseif( relhum(i,k) >= 1.0_r8 ) then
          rhdfda(i,k) = 0.0_r8
      else
         ! Under certain circumstances, rh+ cause cld not to changed
         ! when at an upper limit, or w/ strong subsidence
         if( ( cld2(i,k) - cld(i,k) ) < 1.e-4_r8 ) then
               rhdfda(i,k) = 0.01_r8*relhum(i,k)*1.e+4_r8  
         else
               rhdfda(i,k) = 0.01_r8*relhum(i,k)/(cld2(i,k)-cld(i,k))
         endif
      endif
   enddo
   enddo

 
   ! ---------------------------------------------- !
   ! Stratiform Cloud Macrophysics and Microphysics !
   ! ---------------------------------------------- !

   call t_startf('microp_driver_microphys')

   lchnk  = state1%lchnk
   ncol   = state1%ncol
   rdtime = 1._r8/dtime

 ! Define fractional amount of stratus condensate and precipitation in ice phase.
 ! This uses a ramp ( -30 ~ -10 for fice, -5 ~ 0 for fsnow ). 
 ! The ramp within convective cloud may be different

   call cldwat_fice( ncol, state1%t, fice, fsnow )

   ptend_loc%name         = 'cldwat'
   ptend_loc%ls           = .true.
   ptend_loc%lq(1)        = .true.
   ptend_loc%lq(ixcldice) = .true.
   ptend_loc%lq(ixcldliq) = .true.

   ptend_loc%lq(ixnumliq) = .true.
   ptend_loc%lq(ixnumice) = .true.


 ! ------------------------------ !
 ! Liquid Microp_Driver Macrophysics !
 ! ------------------------------ !

   call t_startf('mmacro_pcond')

   zeros(:ncol,:pver)  = 0._r8
   qc(:ncol,:pver) = state1%q(:ncol,:pver,ixcldliq)
   qi(:ncol,:pver) = state1%q(:ncol,:pver,ixcldice)
   nc(:ncol,:pver) = state1%q(:ncol,:pver,ixnumliq)
   ni(:ncol,:pver) = state1%q(:ncol,:pver,ixnumice)

 ! In CAM5, 'microphysical forcing' ( CC_... ) and 'the other advective forcings' ( ttend, ... ) 
 ! are separately provided into the prognostic microp_driver macrophysics scheme. This is an
 ! attempt to resolve in-cloud and out-cloud forcings. 

   if( get_nstep() .le. 1 ) then
       tcwat(:ncol,:pver)   = state1%t(:ncol,:pver)
       qcwat(:ncol,:pver)   = state1%q(:ncol,:pver,1)
       lcwat(:ncol,:pver)   = qc(:ncol,:pver) + qi(:ncol,:pver)
       iccwat(:ncol,:pver)  = qi(:ncol,:pver)
       nlwat(:ncol,:pver)   = nc(:ncol,:pver)
       niwat(:ncol,:pver)   = ni(:ncol,:pver)
       ttend(:ncol,:pver)   = 0._r8
       qtend(:ncol,:pver)   = 0._r8
       ltend(:ncol,:pver)   = 0._r8
       itend(:ncol,:pver)   = 0._r8
       nltend(:ncol,:pver)  = 0._r8
       nitend(:ncol,:pver)  = 0._r8
       CC_T(:ncol,:pver)    = 0._r8
       CC_qv(:ncol,:pver)   = 0._r8
       CC_ql(:ncol,:pver)   = 0._r8
       CC_qi(:ncol,:pver)   = 0._r8
       CC_nl(:ncol,:pver)   = 0._r8
       CC_ni(:ncol,:pver)   = 0._r8
       CC_qlst(:ncol,:pver) = 0._r8
   else
       ttend(:ncol,:pver)   = ( state1%t(:ncol,:pver)   -  tcwat(:ncol,:pver)) * rdtime -   CC_T(:ncol,:pver) 
       qtend(:ncol,:pver)   = ( state1%q(:ncol,:pver,1) -  qcwat(:ncol,:pver)) * rdtime -  CC_qv(:ncol,:pver)
       ltend(:ncol,:pver)   = ( qc(:ncol,:pver) + qi(:ncol,:pver)  -  & 
                                     lcwat(:ncol,:pver) ) * rdtime - (CC_ql(:ncol,:pver) + CC_qi(:ncol,:pver))
       itend(:ncol,:pver)   = ( qi(:ncol,:pver)         - iccwat(:ncol,:pver)) * rdtime -  CC_qi(:ncol,:pver)
       nltend(:ncol,:pver)  = ( nc(:ncol,:pver)         -  nlwat(:ncol,:pver)) * rdtime -  CC_nl(:ncol,:pver)
       nitend(:ncol,:pver)  = ( ni(:ncol,:pver)         -  niwat(:ncol,:pver)) * rdtime -  CC_ni(:ncol,:pver)
   endif
   lmitend(:ncol,:pver) = ltend(:ncol,:pver) - itend(:ncol,:pver)

   t_inout(:ncol,:pver)  =  tcwat(:ncol,:pver) 
   qv_inout(:ncol,:pver) =  qcwat(:ncol,:pver)
   ql_inout(:ncol,:pver) =  lcwat(:ncol,:pver) - iccwat(:ncol,:pver)
   qi_inout(:ncol,:pver) = iccwat(:ncol,:pver)
   nl_inout(:ncol,:pver) =  nlwat(:ncol,:pver)
   ni_inout(:ncol,:pver) =  niwat(:ncol,:pver)

 ! Liquid Microp_Driver Macrophysics.
 ! The main roles of this subroutines are
 ! (1) compute net condensation rate of strayiform liquid ( cmeliq )
 ! (2) compute liquid stratus and ice stratus fractions. 
 ! Note 'ttend...' are advective tendencies except microphysical process while
 !      'CC...'    are microphysical tendencies. 

   call mmacro_pcond( lchnk, ncol, dtime, state1%pmid, state1%pdel,              &
                      t_inout, qv_inout, ql_inout, qi_inout, nl_inout, ni_inout, &                  
                      ttend, qtend, lmitend, itend, nltend, nitend,              &
                      CC_T, CC_qv, CC_ql, CC_qi, CC_nl, CC_ni, CC_qlst,          & 
                      dlf_T, dlf_qv, dlf_ql, dlf_qi, dlf_nl, dlf_ni,             &
                      concld_old, concld, landfrac, snowh,                       &
                      tlat, qvlat, qcten, qiten, ncten, niten,                   &
                      cmeliq, qvadj, qladj, qiadj, qllim, qilim,                 &
                      cld, alst, aist, qlst, qist ) 

 ! Copy of concld/fice to put in physics buffer
 ! Below are used only for convective cloud.

   concld_ql(:ncol,:pver) = concld(:ncol,:pver)
   fice_ql(:ncol,:pver)   = fice(:ncol,:pver)


 ! Compute net stratus fraction using maximum over-lapping assumption

   do k = 1, pver
   do i = 1, ncol
      ast(i,k) = max( alst(i,k), aist(i,k) )      
   enddo
   enddo

   call t_stopf('mmacro_pcond')

   do k = 1, pver
   do i = 1, ncol
      ptend_loc%s(i,k)          =  tlat(i,k)
      ptend_loc%q(i,k,1)        = qvlat(i,k)
      ptend_loc%q(i,k,ixcldliq) = qcten(i,k)
      ptend_loc%q(i,k,ixcldice) = qiten(i,k)
      ptend_loc%q(i,k,ixnumliq) = ncten(i,k)
      ptend_loc%q(i,k,ixnumice) = niten(i,k)
   end do
   end do   

   call outfld( 'MACPDT   ', tlat ,  pcols, lchnk )
   call outfld( 'MACPDQ   ', qvlat,  pcols, lchnk )
   call outfld( 'MACPDLIQ ', qcten,  pcols, lchnk )
   call outfld( 'MACPDICE ', qiten,  pcols, lchnk )
   call outfld( 'CLDVAPADJ', qvadj,  pcols, lchnk )
   call outfld( 'CLDLIQADJ', qladj,  pcols, lchnk )
   call outfld( 'CLDICEADJ', qiadj,  pcols, lchnk )
   call outfld( 'CLDLIQDET', dlf_ql, pcols, lchnk )
   call outfld( 'CLDICEDET', dlf_qi, pcols, lchnk )
   call outfld( 'CLDLIQLIM', qllim,  pcols, lchnk )
   call outfld( 'CLDICELIM', qilim,  pcols, lchnk )


   call outfld( 'CONCLD  ', concld, pcols, lchnk )
   call outfld( 'CLDST   ', cldst,  pcols, lchnk )
   call outfld( 'CNVCLD  ', clc,    pcols, lchnk )

! calculations and outfld calls for CLDLIQSTR, CLDICESTR, CLDLIQCON, CLDICECON for CFMIP

   ! initialize local variables
   mr_ccliq(1:ncol,1:pver)=0._r8
   mr_ccice(1:ncol,1:pver)=0._r8
   mr_lsliq(1:ncol,1:pver)=0._r8
   mr_lsice(1:ncol,1:pver)=0._r8

   do k=1,pver
   do i=1,ncol
      if (cld(i,k) .gt. 0._r8) then
    	mr_ccliq(i,k) = 0._r8   !! not seen by radiation, so setting to 0 
	mr_ccice(i,k) = 0._r8   !! not seen by radiation, so setting to 0
	mr_lsliq(i,k) = state%q(i,k,ixcldliq)
	mr_lsice(i,k) = state%q(i,k,ixcldice)
      else
	mr_ccliq(i,k) = 0._r8
	mr_ccice(i,k) = 0._r8
	mr_lsliq(i,k) = 0._r8
	mr_lsice(i,k) = 0._r8
      end if
   end do
   end do


   call outfld( 'CLDLIQSTR  ', mr_lsliq,    pcols, lchnk )
   call outfld( 'CLDICESTR  ', mr_lsice,    pcols, lchnk )
   call outfld( 'CLDLIQCON  ', mr_ccliq,    pcols, lchnk )
   call outfld( 'CLDICECON  ', mr_ccice,    pcols, lchnk )

 ! Here 'state_eq' is the equlibrium state after macrophysics for potential
 ! use in the radiation scheme later.

   call physics_ptend_sum( ptend_loc, ptend_all, state )
   call physics_update( state1, tend, ptend_loc, dtime )
   call physics_state_copy( state1, state_eq )

  end subroutine macrop_driver_tend


  end module macrop_driver
