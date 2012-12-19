  module microp_driver

  !-------------------------------------------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the CAM interface to the prognostic cloud microphysics
  !
  ! Author: Andrew Gettelman, Cheryl Craig, October 2010
  ! Origin: modified from stratiform.F90 
  !         (Boville 2002, Coleman 2004, Park 2009, Kay 2010)
  !-------------------------------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: gravit  
  use phys_control,  only: phys_getopts

  use perf_mod
  use cam_logfile,   only: iulog

  implicit none
  private
  save

  public :: microp_driver_register, microp_driver_init_cnst, microp_driver_implements_cnst
  public :: microp_driver_init
  public :: microp_driver_tend

  ! ------------------------- !
  ! Private Module Parameters !
  ! ------------------------- !

  ! Choose either 'intermediate' ('inter') or complete ('compl') cloud microphysics 
  ! inter : Microphysics assumes 'liquid stratus frac = ice stratus frac = max( liquid stratus frac, ice stratus frac )'.
  ! compl : Microphysics explicitly treats 'liquid stratus frac .ne. ice stratus frac'  
  ! for CAM5, only 'inter' is functional

    character(len=5), private, parameter :: micro_treatment = 'inter' 

    logical, private :: sub_column = .false. ! True = configure microphysics for sub-columns 
                                                     ! False = use in regular mode w/o sub-columns

  ! -------------------------------- !
  ! End of Private Module Parameters !
  ! -------------------------------- !

  integer, parameter :: ncnstmax = 4                    ! Number of constituents
  integer            :: ncnst 	  		        ! Number of constituents (can vary)
  character(len=8), dimension(ncnstmax), parameter &    ! Constituent names
                     :: cnst_names = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE'/)
  logical            :: use_shfrc                       ! Local copy of flag from convect_shallow_use_shfrc
  character(len=16)  :: microp_scheme                   ! Microphysics scheme

  integer :: &
     cldo_idx     ,&! old cld index in physics buffer
#ifdef MODAL_AERO
     dgnumwet_idx ,&
     dgnum_idx    ,&
#endif
     ixcldliq     ,&! cloud liquid amount index
     ixcldice     ,&! cloud ice amount index
     ixnumliq     ,&! cloud liquid number index
     ixnumice     ,&! cloud ice water index
     rel2_idx     ,&! rel2 index in physics buffer
     rei2_idx     ,&! rei2 index in physics buffer
     ls_flxprc_idx, ls_flxsnw_idx, qcwat_idx    , lcwat_idx    ,&
     iccwat_idx   , nlwat_idx    , niwat_idx    , cc_t_idx     ,&
     cc_qv_idx    , cc_ql_idx    , cc_qi_idx    , cc_nl_idx    ,&
     cc_ni_idx    , cc_qlst_idx  , tcwat_idx    , cld_idx      ,&
     ast_idx      , aist_idx     , alst_idx     , concld_idx   ,&
     rhdfda_idx   , rhu00_idx    , kvh_idx      , tke_idx      ,&
     turbtype_idx , smaw_idx     , fice_idx     , qme_idx      ,&
     prain_idx    , nevapr_idx   , wsedl_idx    , rei_idx      ,&
     rel_idx      , rel_fn_idx   , dei_idx      , mu_idx       ,&
     lambdac_idx  , iciwp_idx    , iclwp_idx    , deiconv_idx  ,&
     muconv_idx   , lambdaconv_idx, iciwpst_idx  , iclwpst_idx  ,&
     iciwpconv_idx, iclwpconv_idx, des_idx      , icswp_idx    ,&
     cldfsnow_idx , rate1_cw2pr_st_idx                         ,&
     ls_mrprc_idx,&
     ls_mrsnw_idx,&
     ls_reffrain_idx,&
     ls_reffsnow_idx,&
     cv_reffliq_idx,&
     cv_reffice_idx  


  contains

  ! ===============================================================================

  subroutine microp_driver_register

  !---------------------------------------------------------------------- !
  !                                                                       !
  ! Register the constituents (cloud liquid and cloud ice) and the fields !
  ! in the physics buffer.                                                !
  !                                                                       !
  !---------------------------------------------------------------------- !

    use constituents, only: cnst_add, pcnst
    use physconst,    only: mwdry, cpair
    use phys_buffer,  only: pbuf_times, pbuf_add

  !-----------------------------------------------------------------------

    call phys_getopts(microp_scheme_out=microp_scheme)

	ncnst = 4

  ! Register cloud water and determine index.

    call cnst_add(cnst_names(1), mwdry, cpair, 0._r8, ixcldliq, &
       longname='Grid box averaged cloud liquid amount')
    call cnst_add(cnst_names(2), mwdry, cpair, 0._r8, ixcldice, &
       longname='Grid box averaged cloud ice amount')


       call cnst_add(cnst_names(3), mwdry, cpair, 0._r8, ixnumliq, &
          longname='Grid box averaged cloud liquid number')
       call cnst_add(cnst_names(4), mwdry, cpair, 0._r8, ixnumice, &
          longname='Grid box averaged cloud ice number')


  ! Request physics buffer space for fields that persist across timesteps.

    call pbuf_add('CLDO',    'global',  1, pver, pbuf_times,    cldo_idx)
    call pbuf_add('REL2',    'global',  1, pver, pbuf_times,    rel2_idx)
    call pbuf_add('REI2',    'global',  1, pver, pbuf_times,    rei2_idx)

  ! Physics buffer variables for convective cloud properties.

    call pbuf_add('FICE',       'physpkg', 1, pver, 1, fice_idx)

    call pbuf_add('QME',        'physpkg', 1, pver, 1, qme_idx)
    call pbuf_add('PRAIN' ,     'physpkg', 1, pver, 1, prain_idx)
    call pbuf_add('NEVAPR' ,    'physpkg', 1, pver, 1, nevapr_idx)

    call pbuf_add('WSEDL',      'physpkg', 1, pver, 1, wsedl_idx)

    call pbuf_add('REI',        'physpkg', 1, pver, 1, rei_idx)
    call pbuf_add('REL',        'physpkg', 1, pver, 1, rel_idx)
    call pbuf_add('REL_FN',     'physpkg', 1, pver, 1, rel_fn_idx)          ! REL at fixed number for indirect rad forcing

    call pbuf_add('DEI',        'physpkg', 1, pver, 1, dei_idx)          ! Mitchell ice effective diameter for radiation
    call pbuf_add('MU',         'physpkg', 1, pver, 1, mu_idx)          ! Size distribution shape parameter for radiation
    call pbuf_add('LAMBDAC',    'physpkg', 1, pver, 1, lambdac_idx)          ! Size distribution shape parameter for radiation
    call pbuf_add('ICIWP',      'physpkg', 1, pver, 1, iciwp_idx)          ! In cloud ice water path for radiation
    call pbuf_add('ICLWP',      'physpkg', 1, pver, 1, iclwp_idx)          ! In cloud liquid water path for radiation

    call pbuf_add('DEICONV',    'physpkg', 1, pver, 1, deiconv_idx)          ! Convective ice effective diameter for radiation
    call pbuf_add('MUCONV',     'physpkg', 1, pver, 1, muconv_idx)          ! Convective size distribution shape parameter for radiation
    call pbuf_add('LAMBDACONV', 'physpkg', 1, pver, 1, lambdaconv_idx)          ! Convective size distribution shape parameter for radiation
    call pbuf_add('ICIWPST',    'physpkg', 1, pver, 1, iciwpst_idx)          ! Stratiform only in cloud ice water path for radiation
    call pbuf_add('ICLWPST',    'physpkg', 1, pver, 1, iclwpst_idx)          ! Stratiform in cloud liquid water path for radiation
    call pbuf_add('ICIWPCONV',  'physpkg', 1, pver, 1, iciwpconv_idx)          ! Convective only in cloud ice water path for radiation
    call pbuf_add('ICLWPCONV',  'physpkg', 1, pver, 1, iclwpconv_idx)          ! Convective in cloud liquid water path for radiation

    call pbuf_add('DES',        'physpkg', 1, pver, 1, des_idx)          ! Snow effective diameter for radiation
    call pbuf_add('ICSWP',      'physpkg', 1, pver, 1, icswp_idx)          ! In cloud snow water path for radiation
    call pbuf_add('CLDFSNOW',   'physpkg', 1, pver ,pbuf_times, cldfsnow_idx) ! Cloud fraction for liquid drops + snow

#ifdef MODAL_AERO
    call pbuf_add('RATE1_CW2PR_ST', 'physpkg', 1, pver, 1, rate1_cw2pr_st_idx)   ! rce 2010/05/01
#endif

    call pbuf_add('LS_FLXPRC',  'physpkg', 1, pverp, 1, ls_flxprc_idx)
    call pbuf_add('LS_FLXSNW',  'physpkg', 1, pverp, 1, ls_flxsnw_idx)

!additional pbuf calls for CAM5 mixing ratio inputs to COSP
    call pbuf_add('LS_MRPRC',    'physpkg', 1, pver, 1, ls_mrprc_idx)
    call pbuf_add('LS_MRSNW',    'physpkg', 1, pver, 1, ls_mrsnw_idx)
    call pbuf_add('LS_REFFRAIN',    'physpkg', 1, pver, 1, ls_reffrain_idx)
    call pbuf_add('LS_REFFSNOW',    'physpkg', 1, pver, 1, ls_reffsnow_idx)
    call pbuf_add('CV_REFFLIQ',    'physpkg', 1, pver, 1, cv_reffliq_idx)
    call pbuf_add('CV_REFFICE',    'physpkg', 1, pver, 1, cv_reffice_idx)

  end subroutine microp_driver_register

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  function microp_driver_implements_cnst(name)

  !----------------------------------------------------------------------------- ! 
  !                                                                              !    
  ! Purpose: return true if specified constituent is implemented by this package !
  !                                                                              !
  ! Author: B. Eaton                                                             !
  !                                                                              ! 
  !----------------------------------------------------------------------------- !
     implicit none
  !-----------------------------Arguments---------------------------------
     character(len=*), intent(in) :: name      ! constituent name
     logical :: microp_driver_implements_cnst     ! return value
  !---------------------------Local workspace-----------------------------
     integer :: m
  !-----------------------------------------------------------------------

     microp_driver_implements_cnst = .false.

     do m = 1, ncnst
        if (name == cnst_names(m)) then
           microp_driver_implements_cnst = .true.
           return
        end if
     end do
  end function microp_driver_implements_cnst

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine microp_driver_init_cnst(name, q, gcid)

  !----------------------------------------------------------------------- !
  !                                                                        !
  ! Initialize the cloud water mixing ratios (liquid and ice), if they are !
  ! not read from the initial file                                         ! 
  !                                                                        !
  !----------------------------------------------------------------------- !
    implicit none
  !---------------------------- Arguments ---------------------------------
    character(len=*), intent(in)  :: name     ! constituent name
    real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
    integer,          intent(in)  :: gcid(:)  ! global column id
  !-----------------------------------------------------------------------

    if ( name == 'CLDLIQ' ) then
       q = 0.0_r8
       return
    else if ( name == 'CLDICE' ) then
       q = 0.0_r8
       return
    else if ( name == 'NUMLIQ' ) then
       q = 0.0_r8
       return
    else if ( name == 'NUMICE' ) then
       q = 0.0_r8
       return
    end if

  end subroutine microp_driver_init_cnst

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine microp_driver_init

  !-------------------------------------------- !
  !                                             !
  ! Initialize the cloud water parameterization !
  !                                             ! 
  !-------------------------------------------- !

    use cldwat,          only: inimc
    use microp_aero,     only: ini_microp_aero
    use cldwat2m_micro,  only: ini_micro
    use cldwat2m_macro,  only: ini_macro
    use constituents,    only: cnst_get_ind, cnst_name, cnst_longname, sflxnam, apcnst, bpcnst
    use cam_history,     only: addfld, add_default, phys_decomp
    use physconst,       only: tmelt, rh2o, rhodair
    use phys_buffer,     only: pbuf_get_fld_idx
    use convect_shallow, only: convect_shallow_use_shfrc
    use dycore,          only: dycore_is
    use phys_control,    only: cam_physpkg_is
#ifdef MODAL_AERO
    use ndrop,           only: activate_init
!   use cam_history,     only: fieldname_len
!   use spmd_utils,      only: masterproc
!   use modal_aero_data, only: cnst_name_cw, &
!                              lmassptr_amode, lmassptrcw_amode, &
!                              nspec_amode, ntot_amode, numptr_amode, numptrcw_amode, ntot_amode
#endif

    integer              :: m, mm
!   logical              :: history_aerosol      ! Output the MAM aerosol tendencies
    logical              :: history_microphysics ! Output variables for microphysics diagnostics package
    logical              :: history_budget       ! Output tendencies and state variables for CAM4
                                                 ! temperature, water vapor, cloud ice and cloud
                                                 ! liquid budgets.
    integer              :: history_budget_histfile_num ! output history file number for budget fields
!#ifdef MODAL_AERO
!    integer                        :: l, lphase, lspec
!    character(len=fieldname_len)   :: tmpname
!    character(len=fieldname_len+3) :: fieldname
!    character(128)                 :: long_name
!    character(8)                   :: unit
!#endif

  !-----------------------------------------------------------------------

!   call phys_getopts( history_aerosol_out        = history_aerosol      , &
!                      history_microphysics_out   = history_microphysics , & 
!                      history_budget_out         = history_budget       , &
!                      history_budget_histfile_num_out = history_budget_histfile_num)

    call phys_getopts( history_microphysics_out   = history_microphysics , & 
                       history_budget_out         = history_budget       , &
                       history_budget_histfile_num_out = history_budget_histfile_num)

  ! Initialization routine for cloud macrophysics and microphysics

	call ini_micro
	call ini_macro
        call ini_microp_aero
!#ifdef MODAL_AERO
!        call activate_init
!#endif

  ! Register history variables

    do m = 1, ncnst
       call cnst_get_ind( cnst_names(m), mm )
       call addfld( cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm)                   , phys_decomp )
       call addfld( sflxnam  (mm), 'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp )
       call add_default( cnst_name(mm), 1, ' ' )
       call add_default( sflxnam  (mm), 1, ' ' )
    enddo

    call addfld (apcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' after physics'  , phys_decomp)
    call addfld (apcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' after physics'  , phys_decomp)
    call addfld (bpcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' before physics' , phys_decomp)
    call addfld (bpcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' before physics' , phys_decomp)

    if( history_budget) then
       call add_default (cnst_name(ixcldliq), history_budget_histfile_num, ' ')
       call add_default (cnst_name(ixcldice), history_budget_histfile_num, ' ')
       call add_default (apcnst   (ixcldliq), history_budget_histfile_num, ' ')
       call add_default (apcnst   (ixcldice), history_budget_histfile_num, ' ')
       call add_default (bpcnst   (ixcldliq), history_budget_histfile_num, ' ')
       call add_default (bpcnst   (ixcldice), history_budget_histfile_num, ' ')
    end if

    call addfld ('CME      ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap within the cloud'                      ,phys_decomp)
    call addfld ('CMEICE   ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap of ice within the cloud'               ,phys_decomp)
    call addfld ('CMELIQ   ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap of liq within the cloud'               ,phys_decomp)
    call addfld ('PRODPREC ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of condensate to precip'              ,phys_decomp)
    call addfld ('EVAPPREC ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling precip'                   ,phys_decomp)
    call addfld ('EVAPSNOW ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling snow'                     ,phys_decomp)
    call addfld ('HPROGCLD ', 'W/kg'    , pver, 'A', 'Heating from prognostic clouds'                          ,phys_decomp)
    call addfld ('FICE     ', 'fraction', pver, 'A', 'Fractional ice content within cloud'                     ,phys_decomp)
    call addfld ('ICWMR    ', 'kg/kg   ', pver, 'A', 'Prognostic in-cloud water mixing ratio'                  ,phys_decomp)
    call addfld ('ICIMR    ', 'kg/kg   ', pver, 'A', 'Prognostic in-cloud ice mixing ratio'                    ,phys_decomp)
    call addfld ('ICWMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus water mixing ratio'                ,phys_decomp)
    call addfld ('ICIMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus ice mixing ratio'                  ,phys_decomp)

  ! MG microphysics diagnostics
    call addfld ('QCSEVAP  ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling cloud water'              ,phys_decomp)
    call addfld ('QISEVAP  ', 'kg/kg/s ', pver, 'A', 'Rate of sublimation of falling cloud ice'                ,phys_decomp)
    call addfld ('QVRES    ', 'kg/kg/s ', pver, 'A', 'Rate of residual condensation term'                      ,phys_decomp)
    call addfld ('CMEIOUT  ', 'kg/kg/s ', pver, 'A', 'Rate of deposition/sublimation of cloud ice'             ,phys_decomp)
    call addfld ('VTRMC    ', 'm/s     ', pver, 'A', 'Mass-weighted cloud water fallspeed'                     ,phys_decomp)
    call addfld ('VTRMI    ', 'm/s     ', pver, 'A', 'Mass-weighted cloud ice fallspeed'                       ,phys_decomp)
    call addfld ('QCSEDTEN ', 'kg/kg/s ', pver, 'A', 'Cloud water mixing ratio tendency from sedimentation'    ,phys_decomp)
    call addfld ('QISEDTEN ', 'kg/kg/s ', pver, 'A', 'Cloud ice mixing ratio tendency from sedimentation'      ,phys_decomp)
    call addfld ('PRAO     ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud water by rain'                        ,phys_decomp)
    call addfld ('PRCO     ', 'kg/kg/s ', pver, 'A', 'Autoconversion of cloud water'                           ,phys_decomp)
    call addfld ('MNUCCCO  ', 'kg/kg/s ', pver, 'A', 'Immersion freezing of cloud water'                       ,phys_decomp)
    call addfld ('MNUCCTO  ', 'kg/kg/s ', pver, 'A', 'Contact freezing of cloud water'                         ,phys_decomp)
    call addfld ('MNUCCDO  ', 'kg/kg/s ', pver, 'A', 'Homogeneous and heterogeneous nucleation from vapor'     ,phys_decomp)
    call addfld ('MNUCCDOhet','kg/kg/s ', pver, 'A', 'Heterogeneous nucleation from vapor'                     ,phys_decomp)
    call addfld ('MSACWIO  ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water from rime-splintering'         ,phys_decomp)
    call addfld ('PSACWSO  ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud water by snow'                        ,phys_decomp)
    call addfld ('BERGSO   ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water to snow from bergeron'         ,phys_decomp)
    call addfld ('BERGO    ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water to cloud ice from bergeron'    ,phys_decomp)
    call addfld ('MELTO    ', 'kg/kg/s ', pver, 'A', 'Melting of cloud ice'                                    ,phys_decomp)
    call addfld ('HOMOO    ', 'kg/kg/s ', pver, 'A', 'Homogeneous freezing of cloud water'                     ,phys_decomp)
    call addfld ('QCRESO   ', 'kg/kg/s ', pver, 'A', 'Residual condensation term for cloud water'              ,phys_decomp)
    call addfld ('PRCIO    ', 'kg/kg/s ', pver, 'A', 'Autoconversion of cloud ice'                             ,phys_decomp)
    call addfld ('PRAIO    ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud ice by rain'                          ,phys_decomp)
    call addfld ('QIRESO   ', 'kg/kg/s ', pver, 'A', 'Residual deposition term for cloud ice'                  ,phys_decomp)
    call addfld ('MNUCCRO  ', 'kg/kg/s ', pver, 'A', 'Heterogeneous freezing of rain to snow'                  ,phys_decomp)
    call addfld ('PRACSO   ', 'kg/kg/s ', pver, 'A', 'Accretion of rain by snow'                               ,phys_decomp)
    call addfld ('MELTSDT  ', 'W/kg    ', pver, 'A', 'Latent heating rate due to melting of snow'              ,phys_decomp)
    call addfld ('FRZRDT   ', 'W/kg    ', pver, 'A', 'Latent heating rate due to homogeneous freezing of rain' ,phys_decomp)
  ! Convective cloud water variables.
    call addfld ('ICIMRCU  ', 'kg/kg   ', pver, 'A', 'Convection in-cloud ice mixing ratio '                   ,phys_decomp)
    call addfld ('ICLMRCU  ', 'kg/kg   ', pver, 'A', 'Convection in-cloud liquid mixing ratio '                ,phys_decomp)	
    call addfld ('ICIMRTOT ', 'kg/kg   ', pver, 'A', 'Total in-cloud ice mixing ratio '                        ,phys_decomp)
    call addfld ('ICLMRTOT ', 'kg/kg   ', pver, 'A', 'Total in-cloud liquid mixing ratio '                     ,phys_decomp)
    call addfld ('ICWMRSH  ', 'kg/kg   ', pver, 'A', 'Shallow Convection in-cloud water mixing ratio '         ,phys_decomp)
    call addfld ('ICWMRDP  ', 'kg/kg   ', pver, 'A', 'Deep Convection in-cloud water mixing ratio '            ,phys_decomp)


    call add_default ('FICE    ', 1, ' ')

    call addfld ('SH_CLD   ', 'fraction', pver, 'A', 'Shallow convective cloud cover'                          ,phys_decomp)
    call addfld ('DP_CLD   ', 'fraction', pver, 'A', 'Deep convective cloud cover'                             ,phys_decomp)
	
    call add_default ('CONCLD  ', 1, ' ')

    call addfld ('AST','fraction',pver, 'A','Stratus cloud fraction',phys_decomp)

  ! History variables for CAM5 microphysics
    call addfld ('MPDT     ', 'W/kg    ', pver, 'A', 'Heating tendency - Morrison microphysics'                ,phys_decomp)
    call addfld ('MPDQ     ', 'kg/kg/s ', pver, 'A', 'Q tendency - Morrison microphysics'                      ,phys_decomp)
    call addfld ('MPDLIQ   ', 'kg/kg/s ', pver, 'A', 'CLDLIQ tendency - Morrison microphysics'                 ,phys_decomp)
    call addfld ('MPDICE   ', 'kg/kg/s ', pver, 'A', 'CLDICE tendency - Morrison microphysics'                 ,phys_decomp)
    call addfld ('MPDW2V   ', 'kg/kg/s ', pver, 'A', 'Water <--> Vapor tendency - Morrison microphysics'       ,phys_decomp)
    call addfld ('MPDW2I   ', 'kg/kg/s ', pver, 'A', 'Water <--> Ice tendency - Morrison microphysics'         ,phys_decomp)
    call addfld ('MPDW2P   ', 'kg/kg/s ', pver, 'A', 'Water <--> Precip tendency - Morrison microphysics'      ,phys_decomp)
    call addfld ('MPDI2V   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Vapor tendency - Morrison microphysics'         ,phys_decomp)
    call addfld ('MPDI2W   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Water tendency - Morrison microphysics'         ,phys_decomp)
    call addfld ('MPDI2P   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Precip tendency - Morrison microphysics'        ,phys_decomp)
    call addfld ('LIQCLDF  ', 'fraction', pver, 'A', 'Stratus Liquid cloud fraction'                           ,phys_decomp)
    call addfld ('ICECLDF  ', 'fraction', pver, 'A', 'Stratus ICE cloud fraction'                              ,phys_decomp)
    call addfld ('IWC      ', 'kg/m3   ', pver, 'A', 'Grid box average ice water content'                      ,phys_decomp)
    call addfld ('LWC      ', 'kg/m3   ', pver, 'A', 'Grid box average liquid water content'                   ,phys_decomp)
    call addfld ('ICWNC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud water number conc'                   ,phys_decomp)
    call addfld ('ICINC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud ice number conc'                     ,phys_decomp)
    call addfld ('EFFLIQ   ', 'Micron  ', pver, 'A', 'Prognostic droplet effective radius'                     ,phys_decomp)
    call addfld ('EFFLIQ_IND','Micron  ', pver, 'A', 'Prognostic droplet effective radius (indirect effect)'   ,phys_decomp)
    call addfld ('EFFICE   ', 'Micron  ', pver, 'A', 'Prognostic ice effective radius'                         ,phys_decomp)
    call addfld ('WSUB     ', 'm/s     ', pver, 'A', 'Diagnostic sub-grid vertical velocity'                   ,phys_decomp)
    call addfld ('WSUBI    ', 'm/s     ', pver, 'A', 'Diagnostic sub-grid vertical velocity for ice'           ,phys_decomp)
    call addfld ('CDNUMC   ', '#/m2    ', 1,    'A', 'Vertically-integrated droplet concentration'             ,phys_decomp)

    if ( history_budget ) then

       call add_default ('EVAPSNOW ', history_budget_histfile_num, ' ')
       call add_default ('EVAPPREC ', history_budget_histfile_num, ' ')
       call add_default ('CMELIQ   ', history_budget_histfile_num, ' ')


          call add_default ('QVRES    ', history_budget_histfile_num, ' ')
          call add_default ('QISEVAP  ', history_budget_histfile_num, ' ')
          call add_default ('QCSEVAP  ', history_budget_histfile_num, ' ')
          call add_default ('QISEDTEN ', history_budget_histfile_num, ' ')
          call add_default ('QCSEDTEN ', history_budget_histfile_num, ' ')
          call add_default ('QIRESO   ', history_budget_histfile_num, ' ')
          call add_default ('QCRESO   ', history_budget_histfile_num, ' ')
          call add_default ('PSACWSO  ', history_budget_histfile_num, ' ')
          call add_default ('PRCO     ', history_budget_histfile_num, ' ')
          call add_default ('PRCIO    ', history_budget_histfile_num, ' ')
          call add_default ('PRAO     ', history_budget_histfile_num, ' ')
          call add_default ('PRAIO    ', history_budget_histfile_num, ' ')
          call add_default ('PRACSO   ', history_budget_histfile_num, ' ')
          call add_default ('MSACWIO  ', history_budget_histfile_num, ' ')
          call add_default ('MPDW2V   ', history_budget_histfile_num, ' ')
          call add_default ('MPDW2P   ', history_budget_histfile_num, ' ')
          call add_default ('MPDW2I   ', history_budget_histfile_num, ' ')
          call add_default ('MPDT     ', history_budget_histfile_num, ' ')
          call add_default ('MPDQ     ', history_budget_histfile_num, ' ')
          call add_default ('MPDLIQ   ', history_budget_histfile_num, ' ')
          call add_default ('MPDICE   ', history_budget_histfile_num, ' ')
          call add_default ('MPDI2W   ', history_budget_histfile_num, ' ')
          call add_default ('MPDI2V   ', history_budget_histfile_num, ' ')
          call add_default ('MPDI2P   ', history_budget_histfile_num, ' ')
          call add_default ('MNUCCTO  ', history_budget_histfile_num, ' ')
          call add_default ('MNUCCRO  ', history_budget_histfile_num, ' ')
          call add_default ('MNUCCCO  ', history_budget_histfile_num, ' ')
          call add_default ('MELTSDT  ', history_budget_histfile_num, ' ')
          call add_default ('MELTO    ', history_budget_histfile_num, ' ')
          call add_default ('HOMOO    ', history_budget_histfile_num, ' ')
          call add_default ('FRZRDT   ', history_budget_histfile_num, ' ')
          call add_default ('CMEIOUT  ', history_budget_histfile_num, ' ')
          call add_default ('BERGSO   ', history_budget_histfile_num, ' ')
          call add_default ('BERGO    ', history_budget_histfile_num, ' ')

    end if

  ! Averaging for cloud particle number and size
    call addfld ('AWNC     ', 'm-3     ', pver, 'A', 'Average cloud water number conc'                         ,phys_decomp)
    call addfld ('AWNI     ', 'm-3     ', pver, 'A', 'Average cloud ice number conc'                           ,phys_decomp)
    call addfld ('AREL     ', 'Micron  ', pver, 'A', 'Average droplet effective radius'                        ,phys_decomp)
    call addfld ('AREI     ', 'Micron  ', pver, 'A', 'Average ice effective radius'                            ,phys_decomp)
  ! Frequency arrays for above
    call addfld ('FREQL    ', 'fraction', pver, 'A', 'Fractional occurance of liquid'                          ,phys_decomp)
    call addfld ('FREQI    ', 'fraction', pver, 'A', 'Fractional occurance of ice'                             ,phys_decomp)

    if( history_microphysics) then
        call add_default ('CDNUMC   ', 1, ' ')
        call add_default ('IWC      ', 1, ' ')
        call add_default ('WSUB     ', 1, ' ')
        call add_default ('FREQL    ', 1, ' ')
        call add_default ('FREQI    ', 1, ' ')
        call add_default ('AREI     ', 1, ' ')
        call add_default ('AREL     ', 1, ' ')
        call add_default ('AWNC     ', 1, ' ')
        call add_default ('AWNI     ', 1, ' ')
    endif

  ! Average cloud top particle size and number (liq, ice) and frequency
    call addfld ('ACTREL   ', 'Micron  ', 1,    'A', 'Average Cloud Top droplet effective radius'              ,phys_decomp)
    call addfld ('ACTREI   ', 'Micron  ', 1,    'A', 'Average Cloud Top ice effective radius'                  ,phys_decomp)
    call addfld ('ACTNL    ', 'Micron  ', 1,    'A', 'Average Cloud Top droplet number'                        ,phys_decomp)
    call addfld ('ACTNI    ', 'Micron  ', 1,    'A', 'Average Cloud Top ice number'                            ,phys_decomp)

    call addfld ('FCTL     ', 'fraction', 1,    'A', 'Fractional occurance of cloud top liquid'                ,phys_decomp)
    call addfld ('FCTI     ', 'fraction', 1,    'A', 'Fractional occurance of cloud top ice'                   ,phys_decomp)

    call add_default ('ICWMR', 1, ' ')
    call add_default ('ICIMR', 1, ' ')

    call addfld ('LS_FLXPRC', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface rain+snow flux', phys_decomp)
    call addfld ('LS_FLXSNW', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface snow flux', phys_decomp)

    call addfld ('REL', 'micron', pver, 'A', 'MG REL stratiform cloud effective radius liquid', phys_decomp)
    call addfld ('REI', 'micron', pver, 'A', 'MG REI stratiform cloud effective radius ice', phys_decomp)
    call addfld ('LS_REFFRAIN', 'micron', pver, 'A', 'ls stratiform rain effective radius', phys_decomp)
    call addfld ('LS_REFFSNOW', 'micron', pver, 'A', 'ls stratiform snow effective radius', phys_decomp)
    call addfld ('CV_REFFLIQ', 'micron', pver, 'A', 'convective cloud liq effective radius', phys_decomp)
    call addfld ('CV_REFFICE', 'micron', pver, 'A', 'convective cloud ice effective radius', phys_decomp)

!#ifdef MODAL_AERO
!! Add dropmixnuc tendencies for all modal aerosol species
!    do m = 1, ntot_amode
!    do lphase = 1, 2
!    do lspec = 0, nspec_amode(m)+1   ! loop over number + chem constituents + water
!       unit = 'kg/m2/s'
!       if (lspec == 0) then   ! number
!          unit = '#/m2/s'
!          if (lphase == 1) then
!             l = numptr_amode(m)
!          else
!             l = numptrcw_amode(m)
!          endif
!       else if (lspec <= nspec_amode(m)) then   ! non-water mass
!          if (lphase == 1) then
!             l = lmassptr_amode(lspec,m)
!          else
!             l = lmassptrcw_amode(lspec,m)
!          endif
!       else   ! water mass
!          cycle
!       end if
!       if (lphase == 1) then
!          tmpname = cnst_name(l)
!       else
!          tmpname = cnst_name_cw(l)
!       end if
!
!       fieldname = trim(tmpname) // '_mixnuc1'
!       long_name = trim(tmpname) // ' dropmixnuc mixnuc column tendency'
!       call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
!       if ( history_aerosol ) then 
!          call add_default( fieldname, 1, ' ' )
!          if ( masterproc ) write(*,'(2a)') 'microp_driver_init addfld - ', fieldname
!       endif
!       
!    end do   ! lspec
!    end do   ! lphase
!    end do   ! m
!#endif

! Retrieve the indices for pbuf entries

    qcwat_idx    = pbuf_get_fld_idx('QCWAT')
    lcwat_idx    = pbuf_get_fld_idx('LCWAT')
    iccwat_idx   = pbuf_get_fld_idx('ICCWAT')
    nlwat_idx    = pbuf_get_fld_idx('NLWAT')
    niwat_idx    = pbuf_get_fld_idx('NIWAT')
    cc_t_idx     = pbuf_get_fld_idx('CC_T')
    cc_qv_idx    = pbuf_get_fld_idx('CC_qv')
    cc_ql_idx    = pbuf_get_fld_idx('CC_ql')
    cc_qi_idx    = pbuf_get_fld_idx('CC_qi')
    cc_nl_idx    = pbuf_get_fld_idx('CC_nl')
    cc_ni_idx    = pbuf_get_fld_idx('CC_ni')
    cc_qlst_idx  = pbuf_get_fld_idx('CC_qlst')
    tcwat_idx    = pbuf_get_fld_idx('TCWAT')
    cld_idx      = pbuf_get_fld_idx('CLD')
    ast_idx      = pbuf_get_fld_idx('AST')
    aist_idx     = pbuf_get_fld_idx('AIST')
    alst_idx     = pbuf_get_fld_idx('ALST')
    concld_idx   = pbuf_get_fld_idx('CONCLD')
    rhdfda_idx   = pbuf_get_fld_idx('RHDFDA')
    rhu00_idx    = pbuf_get_fld_idx('RHU00')
#ifdef MODAL_AERO    
    dgnumwet_idx = pbuf_get_fld_idx('DGNUMWET')
    dgnum_idx    = pbuf_get_fld_idx('DGNUM' )
#endif

    kvh_idx      = pbuf_get_fld_idx('kvh')
    tke_idx      = pbuf_get_fld_idx('tke')
    turbtype_idx = pbuf_get_fld_idx('turbtype')
    smaw_idx     = pbuf_get_fld_idx('smaw')


    return
  end subroutine microp_driver_init

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine microp_driver_tend(                             &
             state, ptend_all, dtime, &
#ifdef MODAL_AERO
             cflx,                                        &
#endif
             rliq, &
             prec_str, snow_str, prec_sed,      &
             snow_sed, prec_pcw, snow_pcw, pbuf, state_eq, cmeliq )

  !-------------------------------------------------------- !  
  !                                                         ! 
  ! Purpose:                                                !
  !                                                         !
  ! Interface to sedimentation, detrain, cloud fraction and !
  !        cloud macro - microphysics subroutines           !
  !                                                         ! 
  ! Author: D.B. Coleman                                    !
  ! Date: Apr 2004                                          !
  !                                                         !
  !-------------------------------------------------------- !

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use ppgrid
  use physics_types,    only: physics_state, physics_ptend, physics_tend
  use physics_types,    only: physics_ptend_init, physics_update, physics_tend_init
  use physics_types,    only: physics_ptend_sum,  physics_state_copy
  use cam_history,      only: outfld
  use cam_history_support, only : fillvalue
  use phys_buffer,      only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx
  use constituents,     only: cnst_get_ind, pcnst
  use cldwat2m_micro,   only: mmicro_pcond
  use microp_aero,      only: microp_aero_ts 
  use physconst,        only: cpair
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr
  use time_manager,     only: is_first_step, get_nstep
  use conv_water,       only: conv_water_4rad
  use abortutils,       only: endrun
#ifdef MODAL_AERO
  use modal_aero_data
#endif
  ! Debug
    use phys_debug_util,  only: phys_debug_col
  ! Debug

  implicit none

  ! Debug
    integer icol
  ! Debug

!
! Parameters
!
  real(r8) pnot                  ! Reference pressure
  parameter (pnot = 1.e5_r8)

  !
  ! Input arguments
  !

  type(physics_state), intent(in)    :: state       ! State variables
  type(physics_ptend), intent(inout)   :: ptend_all   ! Package tendencies
  type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf

  real(r8), intent(in)  :: dtime                    ! Timestep
#ifdef MODAL_AERO
  real(r8), intent(in)  :: cflx(pcols,pcnst)        ! Constituent flux from surface
#endif

  real(r8), intent(in)  :: rliq(pcols)              ! Vertical integral of liquid not yet in q(ixcldliq)

  real(r8), intent(out) :: prec_str(pcols)          ! [Total] Sfc flux of precip from stratiform [ m/s ] 
  real(r8), intent(out) :: snow_str(pcols)          ! [Total] Sfc flux of snow from stratiform   [ m/s ]
  real(r8), intent(out) :: prec_sed(pcols)          ! Surface flux of total cloud water from sedimentation
  real(r8), intent(out) :: snow_sed(pcols)          ! Surface flux of cloud ice from sedimentation
  real(r8), intent(out) :: prec_pcw(pcols)          ! Sfc flux of precip from microphysics [ m/s ]
  real(r8), intent(out) :: snow_pcw(pcols)          ! Sfc flux of snow from microphysics [ m/s ]

  ! Equilibrium state variables at the end of macrophysics
  ! Below 'state_eq' is for future use as the input of radiation'PBL scheme

  type(physics_state), intent(in) :: state_eq   ! Equilibrium state variables at the end of macrophysics

  ! used from macrophysics
  real(r8), intent(in) :: cmeliq(pcols,pver)                      ! Rate of cond-evap of liq within the cloud

  !
  ! Local variables
  !

  type(physics_state)   :: state1                   ! Local copy of the state variable
  type(physics_tend )   :: tend                     ! Physics tendencies (empty, needed for physics_update call)
  type(physics_ptend)   :: ptend_loc                ! Package tendencies

  integer i,k,m
  integer :: lchnk                                  ! Chunk identifier
  integer :: ncol                                   ! Number of atmospheric columns
  integer :: conv_water_in_rad

  ! Physics buffer fields

  integer itim, ifld
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
  real(r8), pointer, dimension(:,:) :: prain        ! Total precipitation (rain + snow)
  real(r8), pointer, dimension(:,:) :: nevapr       ! Evaporation of total precipitation (rain + snow)
  real(r8), pointer, dimension(:,:) :: rel          ! Liquid effective drop radius (microns)
  real(r8), pointer, dimension(:,:) :: rei          ! Ice effective drop size (microns)
  real(r8), pointer, dimension(:,:) :: rel2         ! Liquid effective drop radius (microns)
  real(r8), pointer, dimension(:,:) :: rei2         ! Ice effective drop size (microns)
  real(r8), pointer, dimension(:,:) :: cldo         ! Old cloud fraction
  real(r8), pointer, dimension(:,:) :: kkvh         ! Vertical eddy diffusivity
  real(r8), pointer, dimension(:,:) :: wsedl        ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
#ifdef MODAL_AERO
  real(r8), pointer, dimension(:,:,:) :: dgnumwet   ! Number mode diameter
  real(r8), pointer, dimension(:,:,:) :: dgnum      ! Number mode diameter
  real(r8), pointer, dimension(:,:) :: rate1ord_cw2pr_st   ! 1st order rate for direct conversion of
                                                    ! strat. cloud water to precip (1/s)    ! rce 2010/05/01
#endif
  real(r8) :: shfrc(pcols,pver)                     ! Cloud fraction from shallow convection scheme
  real(r8), pointer, dimension(:,:) :: rel_fn       ! Ice effective drop size at fixed number (indirect effect) (microns)

  real(r8)  rate1cld(pcols,pver)                    ! array to hold rate1ord_cw2pr_st from microphysics

  ! physics buffer fields for COSP simulator
  real(r8), pointer, dimension(:,:) :: mgflxprc     ! MG grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
  real(r8), pointer, dimension(:,:) :: mgflxsnw     ! MG grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)
  real(r8), pointer, dimension(:,:) :: mgmrprc      ! MG grid-box mean mixingratio_large_scale_cloud_rain+snow at interfaces (kg/kg)
  real(r8), pointer, dimension(:,:) :: mgmrsnw      ! MG grid-box mean mixingratio_large_scale_cloud_snow at interfaces (kg/kg)
  real(r8), pointer, dimension(:,:) :: mgreffrain   ! MG diagnostic rain effective radius (um)
  real(r8), pointer, dimension(:,:) :: mgreffsnow   ! MG diagnostic snow effective radius (um)
  real(r8), pointer, dimension(:,:) :: cvreffliq    ! convective cloud liquid effective radius (um)
  real(r8), pointer, dimension(:,:) :: cvreffice    ! convective cloud ice effective radius (um)

  ! physics buffer fields for radiation

  real(r8), pointer, dimension(:,:) :: dei          ! Ice effective diameter (meters) (AG: microns?)
  real(r8), pointer, dimension(:,:) :: mu           ! Size distribution shape parameter for radiation
  real(r8), pointer, dimension(:,:) :: lambdac      ! Size distribution slope parameter for radiation
  real(r8), pointer, dimension(:,:) :: iciwp        ! In-cloud ice water path for radiation
  real(r8), pointer, dimension(:,:) :: iclwp        ! In-cloud liquid water path for radiation
  
  ! For rrtm optics. specificed distribution.

  real(r8) :: mucon                                 ! Convective size distribution shape parameter
  real(r8) :: dcon                                  ! Convective size distribution effective radius (meters)  
  real(r8) :: lamcon                                ! Convective size distribution slope parameter (meters-1)
  real(r8) :: deicon                                ! Convective ice effective diameter (meters)

  ! Physics buffer fields

  real(r8), pointer, dimension(:,:) :: deiconv      ! Ice effective diameter (microns)
  real(r8), pointer, dimension(:,:) :: muconv       ! Size distribution shape parameter for radiation
  real(r8), pointer, dimension(:,:) :: lambdaconv   ! Size distribution slope parameter for radiation
  real(r8), pointer, dimension(:,:) :: iciwpst      ! Stratiform in-cloud ice water path for radiation
  real(r8), pointer, dimension(:,:) :: iclwpst      ! Stratiform in-cloud liquid water path for radiation
  real(r8), pointer, dimension(:,:) :: iciwpconv    ! Convective in-cloud ice water path for radiation
  real(r8), pointer, dimension(:,:) :: iclwpconv    ! Convective in-cloud liquid water path for radiation

  real(r8), pointer, dimension(:,:) :: tke          ! TKE from the moist PBL scheme
  real(r8), pointer, dimension(:,:) :: turbtype     ! Turbulence type from the moist PBL scheme
  real(r8), pointer, dimension(:,:) :: smaw         ! Instability function of momentum from the moist PBL scheme

  ! Convective cloud to the physics buffer for purposes of ql contrib. to radn.

  real(r8), pointer, dimension(:,:) :: fice      ! Cloud ice/water partitioning ratio.

  ! Local variables for in-cloud water quantities adjusted for convective water

  real(r8)  allcld_ice (pcols,pver)                 ! All-cloud cloud ice
  real(r8)  allcld_liq (pcols,pver)                 ! All-cloud liquid

  ! Snow

  real(r8), pointer, dimension(:,:) :: cldfsnow     ! Cloud fraction for liquid+snow
  real(r8), pointer, dimension(:,:) :: icswp        ! In-cloud snow water path
  real(r8), pointer, dimension(:,:) :: des          ! Snow effective diameter (m)
  real(r8)  qsout(pcols,pver)                       ! Snow mixing ratio

  real(r8) :: qrout(pcols,pver)                     ! Rain mixing ratio
  real(r8) :: rflx(pcols,pver+1)                   ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8) :: sflx(pcols,pver+1)                   ! grid-box average snow flux (kg m^-2 s^-1)
  real(r8) :: reff_rain(pcols,pver)                ! rain effective radius (um)
  real(r8) :: reff_snow(pcols,pver)                ! snow effective radius (um)
 
  real(r8)  icecldf(pcols,pver)                     ! Ice cloud fraction
  real(r8)  liqcldf(pcols,pver)                     ! Liquid cloud fraction (combined into cloud)
  real(r8)  icecldf_out(pcols,pver)                 ! Ice cloud fraction
  real(r8)  liqcldf_out(pcols,pver)                 ! Liquid cloud fraction (combined into cloud)

  ! Local variables for microphysics

  real(r8)  rdtime                                  ! 1./dtime
  real(r8)  qtend(pcols,pver)                       ! Moisture tendencies
  real(r8)  ttend(pcols,pver)                       ! Temperature tendencies
  real(r8)  ltend(pcols,pver)                       ! Cloud liquid water tendencies
  real(r8)  evapsnow(pcols,pver)                    ! Local evaporation of snow
  real(r8)  prodsnow(pcols,pver)                    ! Local production of snow
  real(r8)  icimr(pcols,pver)                       ! In cloud ice mixing ratio
  real(r8)  icwmr(pcols,pver)                       ! In cloud water mixing ratio
  real(r8)  icimrst(pcols,pver)                     ! In stratus ice mixing ratio
  real(r8)  icwmrst(pcols,pver)                     ! In stratus water mixing ratio
  real(r8)  icimrst_out(pcols,pver)                 ! In stratus ice mixing ratio
  real(r8)  icwmrst_out(pcols,pver)                 ! In stratus water mixing ratio
  real(r8)  cmeice(pcols,pver)                      ! Rate of cond-evap of ice within the cloud
  real(r8)  temp(pcols)
  real(r8)  res(pcols,pver)

! MG micro diagnostics

  real(r8)  qcsevap(pcols,pver)                     ! Evaporation of falling cloud water
  real(r8)  qisevap(pcols,pver)                     ! Sublimation of falling cloud ice
  real(r8)  qvres(pcols,pver)                       ! Residual condensation term to remove excess saturation
  real(r8)  cmeiout(pcols,pver)                     ! Deposition/sublimation rate of cloud ice
  real(r8)  vtrmc(pcols,pver)                       ! Mass-weighted cloud water fallspeed
  real(r8)  vtrmi(pcols,pver)                       ! Mass-weighted cloud ice fallspeed
  real(r8)  qcsedten(pcols,pver)                    ! Cloud water mixing ratio tendency from sedimentation
  real(r8)  qisedten(pcols,pver)                    ! Cloud ice mixing ratio tendency from sedimentation

  real(r8)  prao(pcols,pver)  
  real(r8)  prco(pcols,pver)  
  real(r8)  mnuccco(pcols,pver)  
  real(r8)  mnuccto(pcols,pver)  
  real(r8)  mnuccdo(pcols,pver)
  real(r8)  mnuccdohet(pcols,pver)
  real(r8)  msacwio(pcols,pver)  
  real(r8)  psacwso(pcols,pver)  
  real(r8)  bergso(pcols,pver)  
  real(r8)  bergo(pcols,pver)  
  real(r8)  melto(pcols,pver)  
  real(r8)  homoo(pcols,pver)  
  real(r8)  qcreso(pcols,pver)  
  real(r8)  prcio(pcols,pver)  
  real(r8)  praio(pcols,pver)  
  real(r8)  qireso(pcols,pver)
  real(r8)  ftem(pcols,pver)
  real(r8)  mnuccro(pcols,pver) 
  real(r8)  pracso (pcols,pver) 
  real(r8)  meltsdt(pcols,pver) 
  real(r8)  frzrdt (pcols,pver) 
  real(r8)  dpdlfliq(pcols,pver)
  real(r8)  dpdlfice(pcols,pver)
  real(r8)  shdlfliq(pcols,pver)
  real(r8)  shdlfice(pcols,pver)
  real(r8)  dpdlft  (pcols,pver)
  real(r8)  shdlft  (pcols,pver)

#ifdef MODAL_AERO
  integer l, lnum, lnumcw, lmass, lmasscw
#endif

  ! Variables for MG microphysics

  real(r8)  dum1,dum2
  real(r8)  qc(pcols,pver)
  real(r8)  qi(pcols,pver)
  real(r8)  nc(pcols,pver)
  real(r8)  ni(pcols,pver)
  real(r8)  icinc(pcols,pver)                       ! In cloud ice number conc
  real(r8)  cdnumc(pcols)                           ! Vertically-integrated droplet concentration
  real(r8)  icwnc(pcols,pver)                       ! In cloud water number conc
  real(r8)  iwc(pcols,pver)                         ! Grid box average ice water content
  real(r8)  lwc(pcols,pver)                         ! Grid box average liquid water content  
  real(r8)  effliq(pcols,pver)                      ! In cloud liq eff rad
  real(r8)  effice(pcols,pver)                      ! In cloud ice eff rad
  real(r8)  effliq_fn(pcols,pver)                   ! In cloud liq eff rad at fixed number concentration	
  real(r8)  wsub(pcols,pver)                        ! Sub-grid vertical velocity (m/s)
  real(r8)  wsubi(pcols,pver)                       ! Sub-grid vertical velocity for ice (m/s)

  ! Output from mmicro_pcond

  real(r8)  tlat(pcols,pver)
  real(r8)  qvlat(pcols,pver)
  real(r8)  qcten(pcols,pver)
  real(r8)  qiten(pcols,pver)
  real(r8)  ncten(pcols,pver)
  real(r8)  niten(pcols,pver)
  real(r8)  effc(pcols,pver)
  real(r8)  effc_fn(pcols,pver)                     ! Liquid effective radius at fixed number (for indirect calc)
  real(r8)  effi(pcols,pver)
  real(r8)  prect(pcols)
  real(r8)  preci(pcols)

  ! Output from microp_aero_ts for aerosol actication
  real(r8)  naai(pcols,pver)      !ice nucleation number
  real(r8)  naai_hom(pcols,pver)  !ice nucleation number (homogeneous)
  real(r8)  npccn(pcols,pver)     !liquid activation number tendency
  real(r8)  rndst(pcols,pver,4)
  real(r8)  nacon(pcols,pver,4)

  ! Averaging arrays for effective radius and number....

  real(r8)  efiout(pcols,pver)
  real(r8)  efcout(pcols,pver)
  real(r8)  ncout(pcols,pver)
  real(r8)  niout(pcols,pver)

  real(r8)  freqi(pcols,pver)
  real(r8)  freql(pcols,pver)

  ! Average cloud top radius & number

  real(r8)  ctrel(pcols)
  real(r8)  ctrei(pcols)
  real(r8)  ctnl(pcols)
  real(r8)  ctni(pcols)
  real(r8)  fcti(pcols)
  real(r8)  fctl(pcols)

  ! Gather mass mixing ratio for all aerosols affecting the climate

  integer               :: naer_all
  real(r8), pointer     :: aermmr1(:,:)
  real(r8), allocatable :: aer_mmr(:,:,:)           ! Aerosol mass mixing ratio

  real(r8)  zeros(pcols,pver)

  real(r8)  alst_mic(pcols,pver)
  real(r8)  aist_mic(pcols,pver)

  real(r8), pointer :: fldcw(:,:)

  ! ======================================================================

  lchnk = state%lchnk
  ncol  = state%ncol

  call phys_getopts( conv_water_in_rad_out = conv_water_in_rad )

  call physics_state_copy(state_eq,state1)             ! Copy state to local state1.
  call physics_ptend_init(ptend_loc)                ! Initialize local ptend type
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

  concld => pbuf(concld_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  cldo => pbuf(cldo_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  rel2 => pbuf(rel2_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  rei2 => pbuf(rei2_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

#ifdef MODAL_AERO

  dgnumwet => pbuf(dgnumwet_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1:ntot_amode)

  dgnum => pbuf(dgnum_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1:ntot_amode)

  rate1ord_cw2pr_st => pbuf(rate1_cw2pr_st_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  if( is_first_step() ) then
     do i=1,pcnst
        fldcw => qqcw_get_field(i,lchnk,.true.)
        if(associated(fldcw)) then
           fldcw(1:pcols,1:pver)          = 1.e-38_r8
        end if
     end do
     dgnumwet(1:pcols,1:pver,1:ntot_amode) = 0.0_r8
     dgnum(1:pcols,1:pver,1:ntot_amode)    = 0.0_r8
     rate1ord_cw2pr_st(1:pcols,1:pver)     = 0.0_r8   
  endif

#endif

! For purposes of convective ql.

  fice => pbuf(fice_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

  qme  => pbuf(qme_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  prain  => pbuf(prain_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  nevapr  => pbuf(nevapr_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  rel  => pbuf(rel_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  rei  => pbuf(rei_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  rel_fn  => pbuf(rel_fn_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  kkvh => pbuf(kvh_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  dei  => pbuf(dei_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  mu  => pbuf(mu_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  lambdac  => pbuf(lambdac_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iciwp  => pbuf(iciwp_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iclwp  => pbuf(iclwp_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  deiconv  => pbuf(deiconv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  muconv  => pbuf(muconv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  lambdaconv  => pbuf(lambdaconv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iciwpst  => pbuf(iciwpst_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iclwpst  => pbuf(iclwpst_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iciwpconv  => pbuf(iciwpconv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iclwpconv  => pbuf(iclwpconv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  des  => pbuf(des_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  icswp  => pbuf(icswp_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  cldfsnow  => pbuf(cldfsnow_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, itim)

  tke => pbuf(tke_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  turbtype => pbuf(turbtype_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  smaw => pbuf(smaw_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  wsedl  => pbuf(wsedl_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

! If first timestep, initialize heatflux....in pbuf at all time levels.

  if( is_first_step() ) then
      kkvh(:,:)     = 0._r8
      tke(:,:)      = 0._r8
      turbtype(:,:) = 0._r8
      smaw(:,:)     = 0._r8
  endif

  ! Assign default size distribution parameters for no-stratiform clouds (convection only)
  ! Also put into physics buffer for possible separate use by radiation

  dcon   = 25.e-6_r8
  mucon  = 5.3_r8
  deicon = 50._r8

  muconv(:,:)     = mucon
  lambdaconv(:,:) = (mucon + 1._r8)/dcon
  deiconv(:,:)    = deicon


   ! ------------------------------------- !
   ! From here, process computation begins ! 
   ! ------------------------------------- !


     ! ----------------------- !
     ! Microp_Driver Microphysics !
     ! ----------------------- !

       ptend_loc%name         = 'microp'
       ptend_loc%ls           = .true.
       ptend_loc%lq(1)        = .true.
       ptend_loc%lq(ixcldliq) = .true.
       ptend_loc%lq(ixcldice) = .true.
       ptend_loc%lq(ixnumliq) = .true.
       ptend_loc%lq(ixnumice) = .true.

#ifndef MODAL_AERO
       call rad_cnst_get_info( 0, naero = naer_all )
       allocate( aer_mmr( pcols, pver, naer_all ) )
       do m = 1, naer_all
          call rad_cnst_get_aer_mmr( 0, m, state1, pbuf, aermmr1 )
          aer_mmr(:ncol,:,m) = aermmr1(:ncol,:)
       enddo
#endif

#ifdef MODAL_AERO
       do m = 1, ntot_amode
          lnum = numptr_amode(m)
          if( lnum > 0 ) then
              ptend_loc%lq(lnum)= .true.
          endif
          do l = 1, nspec_amode(m)
             lmass = lmassptr_amode(l,m)
             ptend_loc%lq(lmass)= .true.
          enddo
       enddo
#endif

      call t_startf('mmicro_pcond')

      zeros(:ncol,:pver)  = 0._r8
      qc(:ncol,:pver) = state1%q(:ncol,:pver,ixcldliq)
      qi(:ncol,:pver) = state1%q(:ncol,:pver,ixcldice)
      nc(:ncol,:pver) = state1%q(:ncol,:pver,ixnumliq)
      ni(:ncol,:pver) = state1%q(:ncol,:pver,ixnumice)

      if( micro_treatment .eq. 'inter' ) then
          alst_mic(:ncol,:pver) = ast(:ncol,:pver)
          aist_mic(:ncol,:pver) = ast(:ncol,:pver)
      elseif( micro_treatment .eq. 'compl' ) then
          alst_mic(:ncol,:pver) = alst(:ncol,:pver)
          aist_mic(:ncol,:pver) = aist(:ncol,:pver)
      endif

      ! calculate aerosol activation (naai for ice, npccn for liquid) and 
      ! dust size (rndst) and number (nacon) for contact nucleation

      call microp_aero_ts ( lchnk, ncol, dtime, state1%t, zeros,                    & 
                         state1%q(1,1,1), qc, qi,                     &
                         nc, ni, state1%pmid, state1%pdel, ast,                     &
                         alst_mic, aist_mic,                                        &
	                 cldo, state1%pint, state1%rpdel, state1%zm, state1%omega,  &
#ifdef MODAL_AERO
                         state1%q, cflx, ptend_loc%q, dgnumwet, dgnum,        &
#else
                         aer_mmr,                                                   &
#endif
                         kkvh, tke, turbtype, smaw, wsub, wsubi,                    &
                         naai,naai_hom, npccn, rndst,nacon)

      ! Call MG Microphysics

      call mmicro_pcond( sub_column, lchnk, ncol, dtime, state1%t,                  &
                         state1%q(1,1,1), qc, qi,                                   &
                         nc, ni, state1%pmid, state1%pdel, ast,                     &
                         alst_mic, aist_mic,                                        &
	                 cldo,                                                      &
                          rate1cld,                                                 & 
                         naai, npccn, rndst,nacon,                                  &
                         tlat, qvlat,                                               &
                         qcten, qiten, ncten, niten, effc,                          &
                         effc_fn, effi, prect, preci,                               & 
                         nevapr, evapsnow,                                          &
                         prain, prodsnow, cmeice, dei, mu,                          &
                         lambdac, qsout, des,                                       &
                         rflx,sflx, qrout,reff_rain,reff_snow,                      &
                         qcsevap, qisevap, qvres, cmeiout,                          &
                         vtrmc, vtrmi, qcsedten, qisedten,                          &
                         prao, prco, mnuccco, mnuccto, msacwio, psacwso,            &
                         bergso, bergo, melto, homoo, qcreso, prcio, praio, qireso, &
                         mnuccro, pracso, meltsdt, frzrdt , mnuccdo )

      do k=1,pver
      do i=1,ncol
          if (naai(i,k) .gt. 0._r8) then
             mnuccdohet(i,k) = mnuccdo(i,k) - (naai_hom(i,k)/naai(i,k))*mnuccdo(i,k)
          else
             mnuccdohet(i,k) = 0._r8
          end if
      end do
      end do

      mgflxprc  => pbuf(ls_flxprc_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,1)
      mgflxsnw  => pbuf(ls_flxsnw_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,1)
      mgflxprc(:ncol,:pverp) = rflx(:ncol,:pverp) + sflx(:ncol,:pverp)
      mgflxsnw(:ncol,:pverp) = sflx(:ncol,:pverp)

      mgmrprc  => pbuf(ls_mrprc_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
      mgmrsnw  => pbuf(ls_mrsnw_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
      mgmrprc(:ncol,:pver) = qrout(:ncol,:pver) + qsout(:ncol,:pver)
      mgmrsnw(:ncol,:pver) = qsout(:ncol,:pver)

      mgreffrain  => pbuf(ls_reffrain_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
      mgreffsnow  => pbuf(ls_reffsnow_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
      mgreffrain(:ncol,:pver) = reff_rain(:ncol,:pver)
      mgreffsnow(:ncol,:pver) = reff_snow(:ncol,:pver)

      cvreffliq  => pbuf(cv_reffliq_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
      cvreffice  => pbuf(cv_reffice_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
      !! initalize as fillvalue
      cvreffliq(:ncol,:pver) = fillvalue
      cvreffice(:ncol,:pver) = fillvalue
      !! calculate effective radius of convective liquid and ice using dcon and deicon (not used by code, not useful for COSP)
      cvreffliq(:ncol,:pver) = 9.0_r8  !! hard-coded as average of hard-coded values used for deep/shallow convective detrainment (near line 1502)
      cvreffice(:ncol,:pver) = 37.0_r8 !! hard-coded as average of hard-coded values used for deep/shallow convective detrainment (near line 1505)

      call outfld( 'LS_REFFRAIN'  , mgreffrain,       pcols, lchnk )
      call outfld( 'LS_REFFSNOW'  , mgreffsnow,       pcols, lchnk )
      call outfld( 'CV_REFFLIQ'   , cvreffliq,       pcols, lchnk )
      call outfld( 'CV_REFFICE'   , cvreffice,       pcols, lchnk )

    ! Reassign rate1 if modal aerosols
#ifdef MODAL_AERO
      rate1ord_cw2pr_st(1:ncol,1:pver)=rate1cld(1:ncol,1:pver)
#endif
    ! Sedimentation velocity for liquid stratus cloud droplet

      wsedl(:ncol,:pver) = vtrmc(:ncol,:pver)

    ! Nominal values for no microp_driver (convective only) cloud.
    ! Convert snow mixing ratio to microns

      do k = 1, pver
      do i = 1, ncol 
         des(i,k) = des(i,k) * 1.e6_r8
         if( ast(i,k) .lt. 1.e-4_r8 ) then
             mu(i,k) = mucon
             lambdac(i,k) = (mucon + 1._r8)/dcon
             dei(i,k) = deicon
         endif
      end do
      end do 

    ! Microphysical tendencies for use in the macrophysics at the next time step

      CC_T(:ncol,:pver)    =  tlat(:ncol,:pver)/cpair
      CC_qv(:ncol,:pver)   = qvlat(:ncol,:pver)
      CC_ql(:ncol,:pver)   = qcten(:ncol,:pver)
      CC_qi(:ncol,:pver)   = qiten(:ncol,:pver)
      CC_nl(:ncol,:pver)   = ncten(:ncol,:pver)
      CC_ni(:ncol,:pver)   = niten(:ncol,:pver)
      CC_qlst(:ncol,:pver) = qcten(:ncol,:pver)/max(0.01_r8,alst_mic(:ncol,:pver))

    ! Net microp_driver condensation rate

      qme(:ncol,:pver) = cmeliq(:ncol,:pver) + cmeiout(:ncol,:pver) 

      call t_stopf('mmicro_pcond')

#ifndef MODAL_AERO
      deallocate(aer_mmr) 
#endif

!   end if


      do k = 1, pver
      do i = 1, ncol
         ptend_loc%s(i,k)          =  tlat(i,k)
         ptend_loc%q(i,k,1)        = qvlat(i,k)
         ptend_loc%q(i,k,ixcldliq) = qcten(i,k)
         ptend_loc%q(i,k,ixcldice) = qiten(i,k)
         ptend_loc%q(i,k,ixnumliq) = ncten(i,k)
         ptend_loc%q(i,k,ixnumice) = niten(i,k)
      enddo
      enddo
    
    ! For precip, accumulate only total precip in prec_pwc and snow_pwc variables.
    ! Other precip output varirables are set to 0
      prec_pcw(:ncol) = prect(:ncol)
      snow_pcw(:ncol) = preci(:ncol)
      prec_sed(:ncol) = 0._r8
      snow_sed(:ncol) = 0._r8
      prec_str(:ncol) = prec_pcw(:ncol) + prec_sed(:ncol) - rliq(:ncol)
      snow_str(:ncol) = snow_pcw(:ncol) + snow_sed(:ncol) 

   ! ------------------------------- !
   ! Update microphysical tendencies !
   ! ------------------------------- !

   call physics_ptend_sum( ptend_loc, ptend_all, state )
   ptend_all%name = 'cldwat'
   call physics_update( state1, tend, ptend_loc, dtime )
   call physics_ptend_init( ptend_loc ) 

   if( micro_treatment .eq. 'inter' ) then
       icecldf(:ncol,:pver) = ast(:ncol,:pver)
       liqcldf(:ncol,:pver) = ast(:ncol,:pver)
   elseif( micro_treatment .eq. 'compl' ) then
       icecldf(:ncol,:pver) = aist(:ncol,:pver)
       liqcldf(:ncol,:pver) = alst(:ncol,:pver)
   endif

   call outfld( 'ICECLDF ', aist,   pcols, lchnk )
   call outfld( 'LIQCLDF ', alst,   pcols, lchnk )
   call outfld( 'AST',      ast,    pcols, lchnk )   


     ! ------------------------------------------------- !
     ! Save equilibrium state variables for macrophysics !        
     ! at the next time step                             !
     ! ------------------------------------------------- !
       do k = 1, pver
          tcwat(:ncol,k)  = state_eq%t(:ncol,k) 
          qcwat(:ncol,k)  = state_eq%q(:ncol,k,1)
          lcwat(:ncol,k)  = state_eq%q(:ncol,k,ixcldliq) + state_eq%q(:ncol,k,ixcldice)
          iccwat(:ncol,k) = state_eq%q(:ncol,k,ixcldice)
          nlwat(:ncol,k)  = state_eq%q(:ncol,k,ixnumliq)
          niwat(:ncol,k)  = state_eq%q(:ncol,k,ixnumice)         
       end do
     ! Effective droplet radius
       rel(:ncol,:pver)        = effc(:ncol,:pver)
       rel_fn(:ncol,:pver)     = effc_fn(:ncol,:pver)	
       rei(:ncol,:pver)        = effi(:ncol,:pver)
       rel2(:ncol,:pver)       = rel(:ncol,:pver) * 0.9071_r8 ! Convect to effective volume radius assuming pgam = 8
       rei2(:ncol,:pver)       = rei(:ncol,:pver) * 0.6057_r8 ! Convect to effective volume radius at pgam = 0 for ice 
     ! ----------------------------------------------------------- ! 
     ! Adjust in-cloud water values to take account of convective  !
     ! in-cloud water. It is used to calculate the values of       !
     ! icwlp and iciwp to pass to the radiation.                   ! 
     ! ----------------------------------------------------------- !
       allcld_ice(:ncol,:pver) = 0._r8 ! Grid-avg all cloud liquid
       allcld_liq(:ncol,:pver) = 0._r8 ! Grid-avg all cloud ice
       if( conv_water_in_rad /= 0 ) then
	   call conv_water_4rad( lchnk, ncol, pbuf, conv_water_in_rad, rei, state1%pdel, &
	                         state1%q(:,:,ixcldliq), state1%q(:,:,ixcldice),         &
                                 allcld_liq, allcld_ice )
	else
	   allcld_liq(:ncol,:) = state1%q(:ncol,:,ixcldliq)  ! Grid-ave all cloud liquid
	   allcld_ice(:ncol,:) = state1%q(:ncol,:,ixcldice)  !           "        ice 
	end if
      ! ------------------------------------------------------------ !
      ! Compute in cloud ice and liquid mixing ratios                !
      ! Note that 'iclwp, iciwp' are used for radiation computation. !
      ! ------------------------------------------------------------ !
        do k = 1, pver
        do i = 1, ncol
            ! Limits for in-cloud mixing ratios consistent with MG microphysics
            ! in-cloud mixing ratio 0.0001 to 0.005 kg/kg
              icimr(i,k)     = min( allcld_ice(i,k) / max(0.0001_r8,cld(i,k)),0.005_r8 )
              icwmr(i,k)     = min( allcld_liq(i,k) / max(0.0001_r8,cld(i,k)),0.005_r8 )
              icimrst(i,k)   = min( state1%q(i,k,ixcldice) / max(0.0001_r8,icecldf(i,k)),0.005_r8 )
              icwmrst(i,k)   = min( state1%q(i,k,ixcldliq) / max(0.0001_r8,liqcldf(i,k)),0.005_r8 )
              icinc(i,k)     = state1%q(i,k,ixnumice) / max(0.0001_r8,icecldf(i,k)) * state1%pmid(i,k) / (287.15_r8*state1%t(i,k))
              icwnc(i,k)     = state1%q(i,k,ixnumliq) / max(0.0001_r8,liqcldf(i,k)) * state1%pmid(i,k) / (287.15_r8*state1%t(i,k))
 	      iwc(i,k)       = allcld_ice(i,k) * state1%pmid(i,k) / (287.15_r8*state1%t(i,k))
	      lwc(i,k)       = allcld_liq(i,k) * state1%pmid(i,k) / (287.15_r8*state1%t(i,k))
              effliq(i,k)    = effc(i,k)
              effliq_fn(i,k) = effc_fn(i,k)
              effice(i,k)    = effi(i,k)
            ! Calculate total cloud water paths in each layer
              iciwp(i,k)     = icimr(i,k) * state1%pdel(i,k) / gravit
              iclwp(i,k)     = icwmr(i,k) * state1%pdel(i,k) / gravit
            ! Calculate microp_driver cloud water paths in each layer
            ! Note: uses stratiform cloud fraction!
              iciwpst(i,k)   = min(state1%q(i,k,ixcldice)/max(0.0001_r8,ast(i,k)),0.005_r8) * state1%pdel(i,k) / gravit
              iclwpst(i,k)   = min(state1%q(i,k,ixcldliq)/max(0.0001_r8,ast(i,k)),0.005_r8) * state1%pdel(i,k) / gravit
            ! Calculate convective in-cloud LWP.
              iclwpconv(i,k) = max(allcld_liq(i,k) - state%q(i,k,ixcldliq),0._r8)/max(0.0001_r8,concld(i,k)) 
              iciwpconv(i,k) = max(allcld_ice(i,k) - state%q(i,k,ixcldice),0._r8)/max(0.0001_r8,concld(i,k)) 
            ! ------------------------------ !
            ! Adjust cloud fraction for snow !
            ! ------------------------------ !
              cldfsnow(i,k) = cld(i,k)
            ! If cloud and only ice ( no convective cloud or ice ), then set to 0.
              if( ( cldfsnow(i,k) .gt. 1.e-4_r8 ) .and. & 
                  ( concld(i,k)   .lt. 1.e-4_r8 ) .and. & 
                  ( state1%q(i,k,ixcldliq) .lt. 1.e-10_r8 ) ) then
                    cldfsnow(i,k) = 0._r8
              endif
            ! If no cloud and snow, then set to 0.25
              if( ( cldfsnow(i,k) .lt. 1.e-4_r8 ) .and. ( qsout(i,k) .gt. 1.e-6_r8 ) ) then 
                    cldfsnow(i,k) = 0.25_r8
              endif
            ! Calculate in-cloud snow water path
              icswp(i,k) = qsout(i,k) / max( 0.0001_r8, cldfsnow(i,k) ) * state1%pdel(i,k) / gravit
        enddo 
        enddo

      ! --------------------- !
      ! History Output Fields !
      ! --------------------- !

      ! Column droplet concentration

        do i = 1, ncol
           cdnumc(i) = 0._r8
           do k = 1, pver
              cdnumc(i) = cdnumc(i) + state1%q(i,k,ixnumliq)*state1%pdel(i,k)/gravit
           end do
        end do

      ! Averaging for new output fields

        efcout(:,:)      = 0._r8
        efiout(:,:)      = 0._r8
        ncout(:,:)       = 0._r8
        niout(:,:)       = 0._r8	
        freql(:,:)       = 0._r8
        freqi(:,:)       = 0._r8
        liqcldf_out(:,:) = 0._r8
        icecldf_out(:,:) = 0._r8
        icwmrst_out(:,:) = 0._r8
        icimrst_out(:,:) = 0._r8
        do k = 1, pver
        do i = 1, ncol
           if( liqcldf(i,k) .gt. 0.01_r8 .and. icwmrst(i,k) .gt. 5.e-5_r8 ) then
               efcout(i,k) = effc(i,k)
               ncout(i,k)  = icwnc(i,k)
               freql(i,k)  = 1._r8
               liqcldf_out(i,k) = liqcldf(i,k)
               icwmrst_out(i,k) = icwmrst(i,k)
           endif
           if( icecldf(i,k) .gt. 0.01_r8 .and. icimrst(i,k) .gt. 1.e-6_r8 ) then
               efiout(i,k) = effi(i,k)
               niout(i,k)  = icinc(i,k)
               freqi(i,k)  = 1._r8
               icecldf_out(i,k) = icecldf(i,k)
               icimrst_out(i,k) = icimrst(i,k)
            endif
        end do
        end do

        call outfld( 'AREL' , efcout,  pcols, lchnk )
        call outfld( 'AREI' , efiout,  pcols, lchnk )
        call outfld( 'AWNC' , ncout,   pcols, lchnk )
        call outfld( 'AWNI' , niout,   pcols, lchnk )
        call outfld( 'FREQL', freql,   pcols, lchnk )
        call outfld( 'FREQI', freqi,   pcols, lchnk )

      ! Cloud top effective radius and number.

        fcti(:)  = 0._r8
        fctl(:)  = 0._r8
        ctrel(:) = 0._r8
        ctrei(:) = 0._r8
        ctnl(:)  = 0._r8
        ctni(:)  = 0._r8
        do i = 1, ncol
        do k = 1, pver
           if( liqcldf(i,k) .gt. 0.01_r8 .and. icwmrst(i,k) .gt. 1.e-7_r8 ) then
               ctrel(i) = effc(i,k)
               ctnl(i)  = icwnc(i,k)
               fctl(i)  = 1._r8
               exit
           endif
           if( icecldf(i,k) .gt. 0.01_r8 .and. icimrst(i,k) .gt. 1.e-7_r8 ) then
               ctrei(i) = effi(i,k)
               ctni(i)  = icinc(i,k)
               fcti(i)  = 1._r8
               exit
           endif
        enddo
        enddo

        call outfld( 'ACTREL'     , ctrel,     pcols, lchnk )
        call outfld( 'ACTREI'     , ctrei,     pcols, lchnk )
        call outfld( 'ACTNL'      , ctnl,      pcols, lchnk )
        call outfld( 'ACTNI'      , ctni,      pcols, lchnk )
        call outfld( 'FCTL'       , fctl,      pcols, lchnk )
        call outfld( 'FCTI'       , fcti,      pcols, lchnk )

        call outfld( 'MPDT'       , tlat,      pcols, lchnk )
        call outfld( 'MPDQ'       , qvlat,     pcols, lchnk )
        call outfld( 'MPDLIQ'     , qcten,     pcols, lchnk )
        call outfld( 'MPDICE'     , qiten,     pcols, lchnk )
        call outfld( 'ICINC'      , icinc,     pcols, lchnk )
        call outfld( 'ICWNC'      , icwnc,     pcols, lchnk )
        call outfld( 'EFFLIQ'     , effliq,    pcols, lchnk )
        call outfld( 'EFFLIQ_IND' , effliq_fn, pcols, lchnk )
        call outfld( 'EFFICE'     , effice,    pcols, lchnk )
        call outfld( 'WSUB'       , wsub,      pcols, lchnk )
        call outfld( 'WSUBI'      , wsubi,     pcols, lchnk )
        call outfld( 'CDNUMC'     , cdnumc,    pcols, lchnk )

        call outfld('LS_FLXPRC', mgflxprc,    pcols, lchnk )
        call outfld('LS_FLXSNW', mgflxsnw,    pcols, lchnk )
        call outfld('REL', rel,    pcols, lchnk )
        call outfld('REI', rei,    pcols, lchnk )

   ! --------------------------------------------- !
   ! General outfield calls for microphysics       !
   ! --------------------------------------------- !

   call outfld( 'IWC'      , iwc,         pcols, lchnk )
   call outfld( 'LWC'      , lwc,         pcols, lchnk )
   call outfld( 'ICIMR'    , icimr,       pcols, lchnk )
   call outfld( 'ICWMR'    , icwmr,       pcols, lchnk )
   call outfld( 'ICIMRST'  , icimrst_out, pcols, lchnk )
   call outfld( 'ICWMRST'  , icwmrst_out, pcols, lchnk )
   call outfld( 'CME'      , qme,         pcols, lchnk )
   call outfld( 'PRODPREC' , prain,       pcols, lchnk )
   call outfld( 'EVAPPREC' , nevapr,      pcols, lchnk )
   call outfld( 'EVAPSNOW' , evapsnow,    pcols, lchnk )
   call outfld( 'QCSEVAP'  , qcsevap,     pcols, lchnk )
   call outfld( 'QISEVAP'  , qisevap,     pcols, lchnk )
   call outfld( 'QVRES'    , qvres,       pcols, lchnk )
   call outfld( 'CMEIOUT'  , cmeiout,     pcols, lchnk )
   call outfld( 'CMELIQ'   , cmeliq,      pcols, lchnk )
   call outfld( 'VTRMC'    , vtrmc,       pcols, lchnk )
   call outfld( 'VTRMI'    , vtrmi,       pcols, lchnk )
   call outfld( 'QCSEDTEN' , qcsedten,    pcols, lchnk )
   call outfld( 'QISEDTEN' , qisedten,    pcols, lchnk )
   call outfld( 'PRAO'     , prao,        pcols, lchnk )
   call outfld( 'PRCO'     , prco,        pcols, lchnk )
   call outfld( 'MNUCCCO'  , mnuccco,     pcols, lchnk )
   call outfld( 'MNUCCTO'  , mnuccto,     pcols, lchnk )
   call outfld( 'MNUCCDO'  , mnuccdo,     pcols, lchnk )
   call outfld( 'MNUCCDOhet', mnuccdohet, pcols, lchnk )
   call outfld( 'MSACWIO'  , msacwio,     pcols, lchnk )
   call outfld( 'PSACWSO'  , psacwso,     pcols, lchnk )
   call outfld( 'BERGSO'   , bergso,      pcols, lchnk )
   call outfld( 'BERGO'    , bergo,       pcols, lchnk )
   call outfld( 'MELTO'    , melto,       pcols, lchnk )
   call outfld( 'HOMOO'    , homoo,       pcols, lchnk )
   call outfld( 'QCRESO'   , qcreso,      pcols, lchnk )
   call outfld( 'PRCIO'    , prcio,       pcols, lchnk )
   call outfld( 'PRAIO'    , praio,       pcols, lchnk )
   call outfld( 'QIRESO'   , qireso,      pcols, lchnk )
   call outfld( 'MNUCCRO'  , mnuccro,     pcols, lchnk )
   call outfld( 'PRACSO'   , pracso ,     pcols, lchnk )
   call outfld( 'MELTSDT'  , meltsdt,     pcols, lchnk )
   call outfld( 'FRZRDT'   , frzrdt ,     pcols, lchnk )

       ftem(:ncol,:pver) =  qcreso(:ncol,:pver)
       call outfld( 'MPDW2V', ftem, pcols, lchnk )
       ftem(:ncol,:pver) =  melto(:ncol,:pver) - mnuccco(:ncol,:pver) - mnuccto(:ncol,:pver) - &
                            bergo(:ncol,:pver) - homoo  (:ncol,:pver) - msacwio(:ncol,:pver)
       call outfld( 'MPDW2I', ftem, pcols, lchnk )
       ftem(:ncol,:pver) = -prao(:ncol,:pver) - prco(:ncol,:pver) - psacwso(:ncol,:pver) - &
                            bergso(:ncol,:pver)
       call outfld( 'MPDW2P', ftem, pcols, lchnk )
       ftem(:ncol,:pver) =  cmeiout(:ncol,:pver) + qireso (:ncol,:pver)
       call outfld( 'MPDI2V', ftem, pcols, lchnk )
       ftem(:ncol,:pver) = -melto(:ncol,:pver) + mnuccco(:ncol,:pver) + mnuccto(:ncol,:pver) + &
                            bergo(:ncol,:pver) + homoo  (:ncol,:pver) + msacwio(:ncol,:pver)
       call outfld( 'MPDI2W', ftem, pcols, lchnk )
       ftem(:ncol,:pver) = -prcio(:ncol,:pver) - praio  (:ncol,:pver)
       call outfld( 'MPDI2P', ftem, pcols, lchnk )

   call t_stopf('microp_driver_microphys')

   end subroutine microp_driver_tend


  end module microp_driver
