
  module mcshallow

  use cam_history,    only:outfld, addfld, phys_decomp
  use error_function, only: erfc
  use cam_logfile,    only: iulog

  implicit none
  private
  save

  public init_mcshallow
  public compute_mcshallow
  public compute_mcshallow_inv
  
  integer , parameter :: r8 = selected_real_kind(12)    !  8 byte real
  real(r8)            :: xlv                            !  latent heat of vaporization
  real(r8)            :: xlf                            !  latent heat of fusion
  real(r8)            :: xls                            !  latent heat of sublimation
  real(r8)            :: cp                             !  specific heat of dry air
  real(r8)            :: zvir                           !  rh2o/rair - 1
  real(r8)            :: r                              !  gas constant for dry air
  real(r8)            :: g                              !  gravitational constant
  real(r8)            :: ep2                            !  mol wgt water vapor / mol wgt dry air 
  real(r8)            :: p00                            !  reference pressure for exner function
  real(r8)            :: rovcp                          !  R/cp

  contains
  
  real(r8) function exnf(pressure)
    real(r8), intent(in)              :: pressure
    exnf = (pressure/p00)**rovcp
    return
  end function exnf


  subroutine init_mcshallow(kind,xlv_in,cp_in,xlf_in,zvir_in,r_in,g_in,ep2_in)

    !------------------------------------------------------------------------- 
    ! Purpose:  
    ! Initialize time independent variables of the shallow convection package.
    !-------------------------------------------------------------------------

    use cam_history,   only:outfld, addfld, phys_decomp
    use ppgrid,        only:pcols, pver, pverp
    implicit none
    integer , intent(in) :: kind       !  kind of reals being passed in
    real(r8), intent(in) :: xlv_in     !  latent heat of vaporization
    real(r8), intent(in) :: xlf_in     !  latent heat of fusion
    real(r8), intent(in) :: cp_in      !  specific heat of dry air
    real(r8), intent(in) :: zvir_in    !  rh2o/rair - 1
    real(r8), intent(in) :: r_in       !  gas constant for dry air
    real(r8), intent(in) :: g_in       !  gravitational constant
    real(r8), intent(in) :: ep2_in     !  mol wgt water vapor / mol wgt dry air 

    ! ------------------------- !
    ! Internal Output Variables !
    ! ------------------------- !

    call addfld('qtflx_Cu','kg/m2/s',pverp,'A','Cumulus qt flux',phys_decomp)
    call addfld('slflx_Cu','J/m2/s',pverp,'A','Cumulus sl flux',phys_decomp)
    call addfld('uflx_Cu','kg/m/s2',pverp,'A','Cumulus u flux',phys_decomp)
    call addfld('vflx_Cu','kg/m/s2',pverp,'A','Cumulus v flux',phys_decomp)

    call addfld('qtten_Cu','kg/kg/s',pver,'A','qt tendency by cumulus convection',phys_decomp)
    call addfld('slten_Cu','J/kg/s',pver,'A','sl tendency by cumulus convection',phys_decomp)
    call addfld('uten_Cu','m/s2',pver,'A','u tendency by cumulus convection',phys_decomp)
    call addfld('vten_Cu','m/s2',pver,'A','v tendency by cumulus convection',phys_decomp)
    call addfld('qvten_Cu','kg/kg/s',pver,'A','qv tendency by cumulus convection',phys_decomp)
    call addfld('qlten_Cu','kg/kg/s',pver,'A','ql tendency by cumulus convection',phys_decomp)
    call addfld('qiten_Cu','kg/kg/s',pver,'A','qi tendency by cumulus convection',phys_decomp)

    call addfld('cbmf_Cu','kg/m2/s',1,'A','Cumulus Base Mass Flux',phys_decomp)
    call addfld('ufrcinvbase_Cu','no',1,'A','Cumulus Fraction at PBL Top',phys_decomp) 
    call addfld('ufrclcl_Cu','no',1,'A','Cumulus Fraction at LCL',phys_decomp)
    call addfld('winvbase_Cu','m/s',1,'A','Cumulus Vertical Velocity at PBL top',phys_decomp)
    call addfld('wlcl_Cu','m/s',1,'A','Cumulus Vertical Velocity at LCL',phys_decomp)
    call addfld('plcl_Cu','Pa',1,'A','LCL of Source Air',phys_decomp)
    call addfld('pinv_Cu','Pa',1,'A','PBL Top Pressure',phys_decomp)
    call addfld('plfc_Cu','Pa',1,'A','LFC of Source Air',phys_decomp)
    call addfld('pbup_Cu','Pa',1,'A','Highest Level of Positive Cu Buoyancy',phys_decomp)
    call addfld('ppen_Cu','Pa',1,'A','Highest Level where Cu W is 0',phys_decomp)
    call addfld('qtsrc_Cu','kg/kg',1,'A','Source Air qt',phys_decomp)
    call addfld('thlsrc_Cu','K',1,'A','Source Air thl',phys_decomp)
    call addfld('thvlsrc_Cu','K',1,'A','Source Air thvl',phys_decomp)
    call addfld('emfkbup_Cu','kg/m2/s',1,'A','Penetrative Mass Flux at kbup',phys_decomp)
    call addfld('cin_Cu','J/kg',1,'A','CIN upto LFC',phys_decomp)
    call addfld('cinlcl_Cu','J/kg',1,'A','CIN upto LCL',phys_decomp)
    call addfld('cbmflimit_Cu','kg/m2/s',1,'A','cbmflimiter',phys_decomp) 
    call addfld('tkeavg_Cu','m2/s2',1,'A','tkeavg_Cu',phys_decomp) 
    call addfld('zinv_Cu','m',1,'A','PBL Top Height',phys_decomp)
    call addfld('rcwp_Cu','kg/m2',1,'A','Cumulus LWP+IWP',phys_decomp)
    call addfld('rlwp_Cu','kg/m2',1,'A','Cumulus LWP',phys_decomp)
    call addfld('riwp_Cu','kg/m2',1,'A','Cumulus IWP',phys_decomp)
    call addfld('tophgt_Cu','m',1,'A','Cumulus Top Height',phys_decomp)

    call addfld('wu_Cu','m/s',pverp,'A','Cumulus Updraft Vertical Velocity',phys_decomp)
    call addfld('ufrc_Cu','no',pverp,'A','Updraft Fractional Area',phys_decomp)
    call addfld('qtu_Cu','kg/kg',pverp,'A','Cumulus Updraft qt',phys_decomp)
    call addfld('thlu_Cu','K',pverp,'A','Cumulus Updraft thl',phys_decomp)
    call addfld('thvu_Cu','K',pverp,'A','Cumulus Updraft thv',phys_decomp)
    call addfld('uu_Cu','m/s',pverp,'A','Cumulus Updraft uwnd',phys_decomp)
    call addfld('vu_Cu','m/s',pverp,'A','Cumulus Updraft vwnd',phys_decomp)
    call addfld('qtu_emf_Cu','kg/kg',pverp,'A','qt of penatratively entrained air',phys_decomp)
    call addfld('thlu_emf_Cu','K',pverp,'A','thl of penatratively entrained air',phys_decomp)
    call addfld('uu_emf_Cu','m/s',pverp,'A','uwnd of penatratively entrained air',phys_decomp)
    call addfld('vu_emf_Cu','m/s',pverp,'A','vwnd of penatratively entrained air',phys_decomp)
    call addfld('umf_Cu','kg/m2/s',pverp,'A','Cumulus Updraft Mass Flux',phys_decomp)
    call addfld('uemf_Cu','kg/m2/s',pverp,'A','Cumulus Net Mass Flux',phys_decomp)
    call addfld('qcu_Cu','kg/kg',pver,'A','Cumulus updraft LWC+IWC',phys_decomp)
    call addfld('qlu_Cu','kg/kg',pver,'A','Cumulus updraft LWC',phys_decomp)
    call addfld('qiu_Cu','kg/kg',pver,'A','Cumulus updraft IWC',phys_decomp)
    call addfld('cufrc_Cu','no',pver,'A','Cumulus cloud fraction',phys_decomp)
    call addfld('fer_Cu','1/m',pver,'A','Cumulus lateral fractional entrainment rate',phys_decomp)
    call addfld('fdr_Cu','1/m',pver,'A','Cumulus lateral fractional detrainment Rate',phys_decomp)

    call addfld('dwten_Cu','kg/kg/s',pver,'A','Expellsion rate of cumulus cloud water to env.',phys_decomp)
    call addfld('diten_Cu','kg/kg/s',pver,'A','Expellsion rate of cumulus ice water to env.',phys_decomp)
    call addfld('qrten_Cu','kg/kg/s',pver,'A','Production rate of rain by cumulus',phys_decomp)
    call addfld('qsten_Cu','kg/kg/s',pver,'A','Production rate of snow by cumulus',phys_decomp)
    call addfld('flxrain_Cu','kg/m2/s',pverp,'A','Rain flux induced by Cumulus',phys_decomp)
    call addfld('flxsnow_Cu','kg/m2/s',pverp,'A','Snow flux induced by Cumulus',phys_decomp)
    call addfld('ntraprd_Cu','kg/kg/s',pver,'A','Net production rate of rain by Cumulus',phys_decomp)
    call addfld('ntsnprd_Cu','kg/kg/s',pver,'A','Net production rate of snow by Cumulus',phys_decomp)

    call addfld('excessu_Cu','no',pver,'A','Updraft Saturation Excess',phys_decomp)
    call addfld('excess0_Cu','no',pver,'A','Environmental Saturation Excess',phys_decomp)
    call addfld('xc_Cu','no',pver,'A','Critical MIxing Ratio',phys_decomp)
    call addfld('aquad_Cu','no',pver,'A','aquad',phys_decomp)
    call addfld('bquad_Cu','no',pver,'A','bquad',phys_decomp)
    call addfld('cquad_Cu','no',pver,'A','cquad',phys_decomp)
    call addfld('bogbot_Cu','no',pver,'A','Cloud Buoyancy at the Bottom Interface',phys_decomp)
    call addfld('bogtop_Cu','no',pver,'A','Cloud Buoyancy at the Top Interface',phys_decomp)

    call addfld('exit_UWCu_Cu','no',1,'A','exit_UWCu',phys_decomp) 
    call addfld('exit_conden_Cu','no',1,'A','exit_conden',phys_decomp) 
    call addfld('exit_klclmkx_Cu','no',1,'A','exit_klclmkx',phys_decomp) 
    call addfld('exit_klfcmkx_Cu','no',1,'A','exit_klfcmkx',phys_decomp) 
    call addfld('exit_ufrc_Cu','no',1,'A','exit_ufrc',phys_decomp) 
    call addfld('exit_wtw_Cu','no',1,'A','exit_wtw',phys_decomp) 
    call addfld('exit_drycore_Cu','no',1,'A','exit_drycore',phys_decomp) 
    call addfld('exit_wu_Cu','no',1,'A','exit_wu',phys_decomp) 
    call addfld('exit_cufilter_Cu','no',1,'A','exit_cufilter',phys_decomp) 
    call addfld('exit_kinv1_Cu','no',1,'A','exit_kinv1',phys_decomp) 
    call addfld('exit_rei_Cu','no',1,'A','exit_rei',phys_decomp) 

    call addfld('limit_shcu_Cu','no',1,'A','limit_shcu',phys_decomp) 
    call addfld('limit_negcon_Cu','no',1,'A','limit_negcon',phys_decomp) 
    call addfld('limit_ufrc_Cu','no',1,'A','limit_ufrc',phys_decomp) 
    call addfld('limit_ppen_Cu','no',1,'A','limit_ppen',phys_decomp) 
    call addfld('limit_emf_Cu','no',1,'A','limit_emf',phys_decomp) 
    call addfld('limit_cinlcl_Cu','no',1,'A','limit_cinlcl',phys_decomp) 
    call addfld('limit_cin_Cu','no',1,'A','limit_cin',phys_decomp) 
    call addfld('limit_cbmf_Cu','no',1,'A','limit_cbmf',phys_decomp) 
    call addfld('limit_rei_Cu','no',1,'A','limit_rei',phys_decomp) 
    call addfld('ind_delcin_Cu','no',1,'A','ind_delcin',phys_decomp) 

    if (kind .ne. r8) then
       write(iulog,*) 'wrong KIND of reals passed to init_mcshallow -- exiting.'
       stop 'init_mcshallow'
    endif
    xlv   = xlv_in
    xlf   = xlf_in
    xls   = xlv + xlf
    cp    = cp_in
    zvir  = zvir_in
    r     = r_in
    g     = g_in
    ep2   = ep2_in
    p00   = 1.e5_r8
    rovcp = r/cp

  end subroutine init_mcshallow

  subroutine compute_mcshallow_inv( mix      , mkx      , iend     , ncnst     , dt     ,            & 
                                    ps0_inv  , zs0_inv  , p0_inv   , z0_inv    , dp0_inv,            &
                                    u0_inv   , v0_inv   , qv0_inv  , ql0_inv   , qi0_inv,            &
                                    t0_inv   , s0_inv   ,                                            &
                                    tke_inv  , cldfrct_inv, concldfrct_inv,  pblh , cush,            & 
                                    umf_inv  , slflx_inv, qtflx_inv,                                 & 
                                    qvten_inv, qlten_inv, qiten_inv,                                 &
                                    sten_inv , uten_inv , vten_inv ,                                 &
                                    qrten_inv, qsten_inv, precip   , snow, cufrc_inv,                &
                                    qcu_inv  , qlu_inv  , qiu_inv  ,                                 &   
                                    cbmf, qc_inv, rliq,  cnt_inv, cnb_inv, qsat, lchnk ) 

    implicit none
    integer , intent(in)    :: lchnk     
    integer , intent(in)    :: mix
    integer , intent(in)    :: mkx
    integer , intent(in)    :: iend
    integer , intent(in)    :: ncnst
    real(r8), intent(in)    :: dt                       !  time step in seconds : 2*delta_t
    real(r8), intent(in)    :: ps0_inv(mix,mkx+1)       !  environmental pressure at full sigma levels
    real(r8), intent(in)    :: zs0_inv(mix,mkx+1)       !  environmental height at full sigma levels
    real(r8), intent(in)    :: p0_inv(mix,mkx)          !  environmental pressure at half sigma levels
    real(r8), intent(in)    :: z0_inv(mix,mkx)          !  environmental height at half sigma levels
    real(r8), intent(in)    :: dp0_inv(mix,mkx)         !  environmental layer pressure thickness
    real(r8), intent(in)    :: u0_inv(mix,mkx)          !  environmental zonal wind
    real(r8), intent(in)    :: v0_inv(mix,mkx)          !  environmental meridional wind
    real(r8), intent(in)    :: qv0_inv(mix,mkx)         !  environmental specific humidity
    real(r8), intent(in)    :: ql0_inv(mix,mkx)         !  environmental liquid water mixing ratio
    real(r8), intent(in)    :: qi0_inv(mix,mkx)         !  environmental ice mixing ratio
    real(r8), intent(in)    :: t0_inv(mix,mkx)          !  environmental temperature
    real(r8), intent(in)    :: s0_inv(mix,mkx)          !  environmental dry static energy
    real(r8), intent(in)    :: tke_inv(mix,mkx+1)       !  turbulent kinetic energy
    real(r8), intent(in)    :: cldfrct_inv(mix,mkx)     !  total cloud fraction of previous time step
    real(r8), intent(in)    :: concldfrct_inv(mix,mkx)  !  total convective ( shallow + deep ) cloud fraction of previous time step
    real(r8), intent(in)    :: pblh(mix)                !  height of PBL
    real(r8), intent(inout) :: cush(mix)                !  convective scale height
    real(r8), intent(out)   :: umf_inv(mix,mkx+1)       !  updraft mass flux at top of layer
    real(r8), intent(out)   :: qvten_inv(mix,mkx)       !  tendency of specific humidity
    real(r8), intent(out)   :: qlten_inv(mix,mkx)       !  tendency of liquid water mixing ratio
    real(r8), intent(out)   :: qiten_inv(mix,mkx)       !  tendency of ice mixing ratio
    real(r8), intent(out)   :: sten_inv(mix,mkx)        !  tendency of static energy
    real(r8), intent(out)   :: uten_inv(mix,mkx)        !  tendency of zonal wind
    real(r8), intent(out)   :: vten_inv(mix,mkx)        !  tendency of meridional wind
    real(r8), intent(out)   :: qrten_inv(mix,mkx)       !  tendency of rain water mixing ratio
    real(r8), intent(out)   :: qsten_inv(mix,mkx)       !  tendency of snow mixing ratio
    real(r8), intent(out)   :: precip(mix)              !  precipitation flux at the surface [m/s]
    real(r8), intent(out)   :: snow(mix)                !  snow flux at the surface [m/s]
    real(r8), intent(out)   :: rliq(mix)                !  vertical integral of qc
    real(r8), intent(out)   :: slflx_inv(mix,mkx+1)     !  updraft liquid static energy flux
    real(r8), intent(out)   :: qtflx_inv(mix,mkx+1)     !  updraft total water flux
    real(r8), intent(out)   :: cufrc_inv(mix,mkx)       !  shallow cumulus cloud fraction
    real(r8), intent(out)   :: qcu_inv(mix,mkx)         !  updraft liquid+ice mixing ratio
    real(r8), intent(out)   :: qlu_inv(mix,mkx)         !  updraft liquid water mixing ratio
    real(r8), intent(out)   :: qiu_inv(mix,mkx)         !  updraft ice mixing ratio
    real(r8), intent(out)   :: qc_inv(mix,mkx)          !
    real(r8), intent(out)   :: cbmf(mix)                !  
    real(r8), intent(out)   :: cnt_inv(mix)             !  
    real(r8), intent(out)   :: cnb_inv(mix)             !  
    integer , external      :: qsat                     !  function pointer to sat vap pressure function


    real(r8)                :: ps0(mix,0:mkx)           !  environmental pressure at full sigma levels
    real(r8)                :: zs0(mix,0:mkx)           !  environmental height at full sigma levels
    real(r8)                :: p0(mix,mkx)              !  environmental pressure at half sigma levels
    real(r8)                :: z0(mix,mkx)              !  environmental height at half sigma levels
    real(r8)                :: dp0(mix,mkx)             !  environmental layer pressure thickness
    real(r8)                :: u0(mix,mkx)              !  environmental zonal wind
    real(r8)                :: v0(mix,mkx)              !  environmental meridional wind
    real(r8)                :: tke(mix,0:mkx)           !  turbulent kinetic energy
    real(r8)                :: cldfrct(mix,mkx)         !  total cloud fraction of previous time step
    real(r8)                :: concldfrct(mix,mkx)      !  total convective cloud fraction previous time step
    real(r8)                :: qv0(mix,mkx)             !  environmental specific humidity
    real(r8)                :: ql0(mix,mkx)             !  environmental liquid water mixing ratio
    real(r8)                :: qi0(mix,mkx)             !  environmental ice mixing ratio
    real(r8)                :: t0(mix,mkx)              !  environmental temperature
    real(r8)                :: s0(mix,mkx)              !  environmental dry static energy
    real(r8)                :: umf(mix,0:mkx)           !  updraft mass flux at top of layer
    real(r8)                :: qvten(mix,mkx)           !  tendency of specific humidity
    real(r8)                :: qlten(mix,mkx)           !  tendency of liquid water mixing ratio
    real(r8)                :: qiten(mix,mkx)           !  tendency of ice mixing ratio
    real(r8)                :: sten(mix,mkx)            !  tendency of static energy
    real(r8)                :: uten(mix,mkx)            !  tendency of zonal wind
    real(r8)                :: vten(mix,mkx)            !  tendency of meridional wind
    real(r8)                :: qrten(mix,mkx)           !  tendency of rain water mixing ratio
    real(r8)                :: qsten(mix,mkx)           !  tendency of snow mixing ratio
    real(r8)                :: slflx(mix,0:mkx)         !  updraft liquid static energy flux
    real(r8)                :: qtflx(mix,0:mkx)         !  updraft total water flux
    real(r8)                :: cufrc(mix,mkx)           !  shallow cumulus cloud fraction
    real(r8)                :: qcu(mix,mkx)             !  updraft condensate water mixing ratio,'qlu(k)+qiu(k)'
    real(r8)                :: qlu(mix,mkx)             !  updraft liquid water mixing ratio
    real(r8)                :: qiu(mix,mkx)             !  updraft ice mixing ratio
    real(r8)                :: qc(mix,mkx)              !  
    real(r8)                :: cnt(mix)                 !  
    real(r8)                :: cnb(mix)                 !   
    integer                 :: k                        !  vertical index for local fields 
    integer                 :: k_inv                    !  vertical index for incoming fields

    do k = 1, mkx
       k_inv               = mkx + 1 - k
       p0(:iend,k)         = p0_inv(:iend,k_inv)
       u0(:iend,k)         = u0_inv(:iend,k_inv)
       v0(:iend,k)         = v0_inv(:iend,k_inv)
       z0(:iend,k)         = z0_inv(:iend,k_inv)
       dp0(:iend,k)        = dp0_inv(:iend,k_inv)
       qv0(:iend,k)        = qv0_inv(:iend,k_inv)
       ql0(:iend,k)        = ql0_inv(:iend,k_inv)
       qi0(:iend,k)        = qi0_inv(:iend,k_inv)
       t0(:iend,k)         = t0_inv(:iend,k_inv)
       s0(:iend,k)         = s0_inv(:iend,k_inv)
       cldfrct(:iend,k)    = cldfrct_inv(:iend,k_inv)
       concldfrct(:iend,k) = concldfrct_inv(:iend,k_inv)
    end do
    
    do k = 0, mkx
       k_inv               = mkx + 1 - k
       ps0(:iend,k)        = ps0_inv(:iend,k_inv)
       zs0(:iend,k)        = zs0_inv(:iend,k_inv)
       tke(:iend,k)        = tke_inv(:iend,k_inv)
    end do

    call compute_mcshallow( mix  , mkx  , iend  , ncnst , dt   ,          &
                            ps0  , zs0  , p0    , z0    , dp0  ,          &
                            u0   , v0   , qv0   , ql0   , qi0  ,          & 
                            t0   , s0   ,                                 & 
                            tke  , cldfrct, concldfrct, pblh, cush  ,     & 
                            umf  , slflx, qtflx ,                         &  
                            qvten, qlten, qiten ,                         & 
                            sten , uten , vten  ,                         &
                            qrten, qsten, precip, snow, cufrc,            &
                            qcu  , qlu  , qiu   ,                         &
                            cbmf, qc, rliq , cnt, cnb, qsat, lchnk )

    ! Reverse cloud top/base interface indices

       cnt_inv(:iend) = mkx + 1 - cnt(:iend)
       cnb_inv(:iend) = mkx + 1 - cnb(:iend)

    do k = 0, mkx
       k_inv                  = mkx + 1 - k
       umf_inv(:iend,k_inv)   = umf(:iend,k)                 !  updraft mass flux at top of layer
       slflx_inv(:iend,k_inv) = slflx(:iend,k)               !  updraft liquid static energy flux
       qtflx_inv(:iend,k_inv) = qtflx(:iend,k)               !  updraft total water flux
    end do

    do k = 1, mkx
       k_inv                         = mkx + 1 - k
       qvten_inv(:iend,k_inv)        = qvten(:iend,k)        ! tendency of specific humidity
       qlten_inv(:iend,k_inv)        = qlten(:iend,k)        ! tendency of liquid water mixing ratio
       qiten_inv(:iend,k_inv)        = qiten(:iend,k)        ! tendency of ice mixing ratio
       sten_inv(:iend,k_inv)         = sten(:iend,k)         ! tendency of static energy
       uten_inv(:iend,k_inv)         = uten(:iend,k)         ! tendency of zonal wind
       vten_inv(:iend,k_inv)         = vten(:iend,k)         ! tendency of meridional wind
       qrten_inv(:iend,k_inv)        = qrten(:iend,k)        ! tendency of rain water mixing ratio
       qsten_inv(:iend,k_inv)        = qsten(:iend,k)        ! tendency of snow mixing ratio
       cufrc_inv(:iend,k_inv)        = cufrc(:iend,k)        ! shallow cumulus cloud fraction
       qcu_inv(:iend,k_inv)          = qcu(:iend,k)          ! updraft liquid+ice water mixing ratio
       qlu_inv(:iend,k_inv)          = qlu(:iend,k)          ! updraft liquid water mixing ratio
       qiu_inv(:iend,k_inv)          = qiu(:iend,k)          ! updraft ice water mixing ratio
       qc_inv(:iend,k_inv)           = qc(:iend,k)           ! reserved "qlten+qiten' associated with cumulus detrained water 
    end do
    
  end subroutine compute_mcshallow_inv


  subroutine compute_mcshallow( mix      , mkx      , iend      , ncnst     , dt     ,           &
                                ps0_in   , zs0_in   , p0_in     , z0_in     , dp0_in ,           &
                                u0_in    , v0_in    , qv0_in    , ql0_in    , qi0_in ,           &
                                t0_in    , s0_in    ,                                            &
                                tke_in   , cldfrct_in, concldfrct_in,  pblh_in  , cush_inout,    & 
                                umf_out  , slflx_out, qtflx_out ,                                &
                                qvten_out, qlten_out, qiten_out ,                                & 
                                sten_out , uten_out , vten_out  ,                                &
                                qrten_out, qsten_out, precip_out, snow_out, cufrc_out,           &
                                qcu_out  , qlu_out  , qiu_out   ,                                &
                                cbmf_out, qc_out, rliq_out, cnt_out, cnb_out, qsat, lchnk )

    !======================================================!
    !                                                      ! 
    !     SHALLOW CONVECTION SCHEME                        !
    !     Described in McCaa, Bretherton, and Grenier:     !
    !     (submitted to MWR, December 2001)                !
    !                                                      !
    !     Modified by Sungsu Park. Oct.2005.               !
    !                                                      !
    !======================================================!
 


    !
    ! Input-Output variables
    !

    use cam_history,   only:outfld, addfld, phys_decomp
    implicit none

    integer , intent(in)    :: lchnk
    integer , intent(in)    :: mix
    integer , intent(in)    :: mkx
    integer , intent(in)    :: iend
    integer , intent(in)    :: ncnst
    real(r8), intent(in)    :: dt                             !  time step in seconds: 2*delta_t
    real(r8), intent(in)    :: ps0_in(mix,0:mkx)              !  environmental pressure at full sigma levels
    real(r8), intent(in)    :: zs0_in(mix,0:mkx)              !  environmental height at full sigma levels
    real(r8), intent(in)    :: p0_in(mix,mkx)                 !  environmental pressure at half sigma levels
    real(r8), intent(in)    :: z0_in(mix,mkx)                 !  environmental height at half sigma levels
    real(r8), intent(in)    :: dp0_in(mix,mkx)                !  environmental layer pressure thickness
    real(r8), intent(in)    :: u0_in(mix,mkx)                 !  environmental zonal wind
    real(r8), intent(in)    :: v0_in(mix,mkx)                 !  environmental meridional wind
    real(r8), intent(in)    :: qv0_in(mix,mkx)                !  environmental specific humidity
    real(r8), intent(in)    :: ql0_in(mix,mkx)                !  environmental liquid water mixing ratio
    real(r8), intent(in)    :: qi0_in(mix,mkx)                !  environmental ice mixing ratio
    real(r8), intent(in)    :: t0_in(mix,mkx)                 !  environmental temperature
    real(r8), intent(in)    :: s0_in(mix,mkx)                 !  environmental dry static energy
    real(r8), intent(in)    :: tke_in(mix,0:mkx)              !  turbulent kinetic energy
    real(r8), intent(in)    :: cldfrct_in(mix,mkx)            !  total cloud fraction at the previous time step
    real(r8), intent(in)    :: concldfrct_in(mix,mkx)         !  total convective cloud fraction at the previous time step
    real(r8), intent(in)    :: pblh_in(mix)                   !  height of PBL
    real(r8), intent(inout) :: cush_inout(mix)                !  convective scale height
    real(r8), intent(out)   :: umf_out(mix,0:mkx)             !  updraft mass flux at top of layer
    real(r8), intent(out)   :: qvten_out(mix,mkx)             !  tendency of specific humidity
    real(r8), intent(out)   :: qlten_out(mix,mkx)             !  tendency of liquid water mixing ratio
    real(r8), intent(out)   :: qiten_out(mix,mkx)             !  tendency of ice mixing ratio
    real(r8), intent(out)   :: sten_out(mix,mkx)              !  tendency of static energy
    real(r8), intent(out)   :: uten_out(mix,mkx)              !  tendency of zonal wind
    real(r8), intent(out)   :: vten_out(mix,mkx)              !  tendency of meridional wind
    real(r8), intent(out)   :: qrten_out(mix,mkx)             !  tendency of rain water mixing ratio
    real(r8), intent(out)   :: qsten_out(mix,mkx)             !  tendency of snow mixing ratio
    real(r8), intent(out)   :: precip_out(mix)                !  precipitation rate at surface [m/s]
    real(r8), intent(out)   :: snow_out(mix)                  !  snow rate at surface [m/s]
    real(r8), intent(out)   :: slflx_out(mix,0:mkx)           !  updraft/pen.entrainment liquid static energy flux
    real(r8), intent(out)   :: qtflx_out(mix,0:mkx)           !  updraft/pen.entrainment total water flux
    real(r8), intent(out)   :: cufrc_out(mix,mkx)             !  shallow cumulus cloud fraction
    real(r8), intent(out)   :: qcu_out(mix,mkx)               !  updraft condensate water mixing ratio, 'qlu(k)+qiu(k)'
    real(r8), intent(out)   :: qlu_out(mix,mkx)               !  updraft liquid water mixing ratio
    real(r8), intent(out)   :: qiu_out(mix,mkx)               !  updraft ice water mixing ratio
    real(r8), intent(out)   :: cbmf_out(mix)                  !  
    real(r8), intent(out)   :: qc_out(mix,mkx)                !  [kg/kg/s]
    real(r8), intent(out)   :: rliq_out(mix)                  !  [m/s]
    real(r8), intent(out)   :: cnt_out(mix)                   !  
    real(r8), intent(out)   :: cnb_out(mix)                   !  
    integer , external      :: qsat 
    real(r8)                   qtten_out(mix,mkx)             !  tendency of qt
    real(r8)                   slten_out(mix,mkx)             !  tendency of sl
    real(r8)                   ufrc_out(mix,0:mkx)            !  updraft fractional area
    real(r8)                   uflx_out(mix,0:mkx)            !  updraft/pen.entrainment zonal momentum flux
    real(r8)                   vflx_out(mix,0:mkx)            !  updraft/pen.entrainment meridional momentum flux
    real(r8)                   fer_out(mix,mkx)               !  lateral entrainment rate
    real(r8)                   fdr_out(mix,mkx)               !  lateral detrainment rate
    real(r8)                   cinh_out(mix)                  !  convective inhibition upto LFC (CIN) ([J/kg])
    
    !
    ! One-dimensional variables at each grid point
    !

    ! 1. Input variables

    real(r8)    ps0(0:mkx)         !  environmental pressure at full sigma levels
    real(r8)    zs0(0:mkx)         !  environmental height at full sigma levels
    real(r8)    p0(mkx)            !  environmental pressure at half sigma levels
    real(r8)    z0(mkx)            !  environmental height at half sigma levels
    real(r8)    dp0(mkx)           !  environmental layer pressure thickness
    real(r8)    u0(mkx)            !  environmental zonal wind
    real(r8)    v0(mkx)            !  environmental meridional wind
    real(r8)    tke(0:mkx)         !  turbulent kinetic energy
    real(r8)    cldfrct(mkx)       !  total cloud fraction at the previous time step
    real(r8)    concldfrct(mkx)    !  total convective cloud fraction at the previous time step
    real(r8)    qv0(mkx)           !  environmental specific humidity
    real(r8)    ql0(mkx)           !  environmental liquid water mixing ratio
    real(r8)    qi0(mkx)           !  environmental ice mixing ratio
    real(r8)    t0(mkx)            !  environmental temperature
    real(r8)    s0(mkx)            !  environmental dry static energy
    real(r8)    pblh               !  height of PBL
    real(r8)    cush               !  convective scale height

    ! 2. Environmental variables directly from the input variables

    real(r8)    qt0(mkx)           !  environmental total water mixing ratio
    real(r8)    thl0(mkx)          !  environmental liquid potential temperature
    real(r8)    thvl0(mkx)         !  environmental liquid virtual potential temperature
    real(r8)    ssqt0(mkx)         !  slope of environmental total water mixing ratio
    real(r8)    ssthl0(mkx)        !  slope of environmental liquid potential temperature
    real(r8)    ssu0(mkx)          !  environmental zonal wind speed vertical gradient
    real(r8)    ssv0(mkx)          !  environmental meridional wind speed vertical gradient
    real(r8)    thv0bot(mkx)       !  environmental virtual potential temperature, bottom of layer
    real(r8)    thv0top(mkx)       !  environmental virtual potential temperature, top of layer
    real(r8)    thvl0bot(mkx)      !  environmental liquid virtual potential temperature, bottom of layer
    real(r8)    thvl0top(mkx)      !  environmental liquid virtual potential temperature, top of layer
    real(r8)    exn0(mkx)          !  exner function at midpoints
    real(r8)    exns0(0:mkx)       !  exner function at interfaces

   ! 3. Variables associated with cumulus convection

    real(r8)    umf(0:mkx)              !  updraft mass flux at top of layer
    real(r8)    emf(0:mkx)              !  updraft mass flux at top of layer
    real(r8)    qvten(mkx)              !  tendency of specific humidity
    real(r8)    qlten(mkx)              !  tendency of liquid water mixing ratio
    real(r8)    qiten(mkx)              !  tendency of ice mixing ratio
    real(r8)    sten(mkx)               !  tendency of static energy
    real(r8)    uten(mkx)               !  tendency of zonal wind
    real(r8)    vten(mkx)               !  tendency of meridional wind
    real(r8)    qrten(mkx)              !  tendency of rain water mixing ratio
    real(r8)    qsten(mkx)              !  tendency of snow mixing ratio
    real(r8)    precip                  !  precipitation rate at the surface [m/s]
    real(r8)    snow                    !  snow rate at the surface [m/s]
    real(r8)    slflx(0:mkx)            !  updraft liquid static energy flux
    real(r8)    qtflx(0:mkx)            !  updraft total water flux
    real(r8)    cufrc(mkx)              !  shallow cumulus cloud fraction
    real(r8)    qcu(mkx)                !  updraft condensate water mixing ratio, 'qlu(k)+qiu(k)'
    real(r8)    qlu(mkx)                !  updraft liquid water mixing ratio
    real(r8)    qiu(mkx)                !  updraft ice water mixing ratio
    real(r8)    uflx(0:mkx)             !  flux of zonal momentum due to convection
    real(r8)    vflx(0:mkx)             !  flux of meridional momentum due to convection
    real(r8)    dwten(mkx)              !  detrained water tendency from cumulus updraft
    real(r8)    diten(mkx)              !  detrained ice   tendency from cumulus updraft 
    real(r8)    fer(mkx)                !  fractional lateral entrainment rate
    real(r8)    fdr(mkx)                !  fractional lateral detrainment rate
    real(r8)    uf(mkx)                 !  zonal wind at end of time step
    real(r8)    vf(mkx)                 !  meridional wind at end of time step
    real(r8)    qc(mkx)                 !  [kg/kg/s] reserved 'qlten+qiten' due to detrained 'cloud water + cloud ice' (without rain-snow contribution)
    real(r8)    qc_l(mkx)
    real(r8)    qc_i(mkx)
    real(r8)    rliq                    !  [m/s] vertical integral of qc 
    real(r8)    cnt, cnb
    real(r8)    qtten(mkx)              !  tendency of qt
    real(r8)    slten(mkx)              !  tendency of sl
    real(r8)    ufrc(0:mkx)             !  updraft fractional area
    
    !----- Variables used for the calculation of condensation sink associated with compensating subsidence

    real(r8)    uemf(0:mkx)             !  net updraft mass flux at the interface of the top of layer ( emf + umf )
    real(r8)    comsub(mkx)             !  compensating subsidence at layer mid-point ( unit of mass flux, umf )
    real(r8)    qcten_sink(mkx)         !  condensate sink ( > 0 -> condensate is generated by subsidence )
    real(r8)    qlten_sink(mkx)         !  liquid condensate sink
    real(r8)    qiten_sink(mkx)         !  ice    condensate sink 
    real(r8)    qs0(mkx)                !  saturation specific humidity at layer mid-point

    !----- Variables describing cumulus updraft

    real(r8)    wu(0:mkx)               !  Updraft vertical velocity at top of layer
    real(r8)    thlu(0:mkx)             !  Updraft liquid potential temperature at top of layer
    real(r8)    uu(0:mkx)               !  Updraft zonal wind speed
    real(r8)    vu(0:mkx)               !  Updraft meridional wind speed
    real(r8)    qtu(0:mkx)              !  Updraft total water at top of layer
    real(r8)    thvu(0:mkx)             !  Updraft virtual potential temperature at top of layer
    real(r8)    rei(mkx)                !  Updraft mixing rate of with environment

    !----- Variables describing conservative scalars of entraining downdrafts  at the 
    !      entraining interfaces, i.e., 'kbup <= k < kpen-1'. At the other interfaces,
    !      belows are simply set to equal to those of updraft for simplicity - but it
    !      does not influence numerical calculation.

    real(r8)    thlu_emf(0:mkx)         !  Downdraft liquid potential temperature at entraining interfaces
    real(r8)    qtu_emf(0:mkx)          !  Downdraft total water at entraining interfaces
    real(r8)    uu_emf(0:mkx)           !  Downdraft zonal wind at entraining interfaces
    real(r8)    vu_emf(0:mkx)           !  Downdraft meridional wind at entraining interfaces
    
    !----- Variables associated with evaporations of 'rain' and 'snow' in bulk microphysics

    real(r8)    flxrain(0:mkx)          !  Downward rain flux at each interface [ kg/m2/s ]
    real(r8)    flxsnow(0:mkx)          !  Downward snow flux at each interface [ kg/m2/s ]
    real(r8)    ntraprd(mkx)            !  Net production rate of rain in each layer [ kg/kg/s ]
    real(r8)    ntsnprd(mkx)            !  Net production rate of snow in each layer [ kg/kg/s ]
    real(r8)    flxsntm                 !  Downward snow flux at the top of each layer after melting [ kg/m2/s ]
    real(r8)    snowmlt                 !  Snow melting tendency [ kg/kg/s ]
    real(r8)    subsat                  !  Sub-saturation ratio (1-qv/qs) [ no unit ]
    real(r8)    evprain                 !  Evaporation rate of rain [ kg/kg/s ]
    real(r8)    evpsnow                 !  Evaporation rate of snow [ kg/kg/s ]
    real(r8)    evplimit_rain           !  Limiter of 'evprain' [ kg/kg/s ]
    real(r8)    evplimit_snow           !  Limiter of 'evpsnow' [ kg/kg/s ]
    real(r8)    evpint_rain             !  Vertically-integrated evaporative flux of rain [ kg/m2/s ]
    real(r8)    evpint_snow             !  Vertically-integrated evaporative flux of snow [ kg/m2/s ]
    real(r8)    kevp                    !  Evaporative efficiency [ complex unit ]

    !----- Other internal variables

    integer     kk
    integer     id_check
    logical     id_exit   
    logical     shallowCu               !  Used for testing whether cumulus updraft can overcome the 
                                        !  buoyancy barrier just above the PBL top.
    real(r8)    thlsrc,qtsrc,usrc,vsrc,uplus,vplus
    real(r8)    plcl,plfc,prel,wrel
    real(r8)    frc_rasn, frc_rasn_remain
    real(r8)    ee2,ud2,wtw,wtwb,wtwh
    real(r8)    cldhgt,dpsum
    real(r8)    xc,xc_2                 !  fraction of environmental air in a neutrally buoyant mixture
    real(r8)    cin                     !  convective inhibition (m2/s2)
    integer     k                       !  release level (half-sigma level just below lcl)
    integer     i
    integer     leff
    integer     kmin                    !  layer where 'thvl0' is minimum within the PBL
    integer     klcl                    !  layer containing LCL of source air
    integer     kinv                    !  inversion layer with PBL top interface as a lower interface
    integer     krel                    !  release layer where buoyancy sorting mixing occurs for the first time
    integer     klfc                    !  LFC layer of cumulus source air
    integer     kbup                    !  top layer in which cloud buoyancy is positive both at the top and bottom interfaces
    integer     kpen                    !  highest layer with positive updraft vertical velocity - top layer cumulus can reach
    integer     iteration
    integer     kp1,km1,m
    real(r8) :: scaleh,nu,tkeavg,thvusrc,qt0bot,qse,cridis,thlxsat,qtxsat,thvxsat,x_cu,x_en,thv_x0,thv_x1
    real(r8) :: thj,qvj, qlj,qij,tscaleh,thlu_top,qtu_top,exntop 
    real(r8) :: cinlcl,rbuoy,rdrag,pe,dpe,dpeh,exne,thvebot,thle,qte,ue,ve,thlue,qtue,wue,thluh,qtuh
    real(r8) :: cbmf,wexp,wcrit,sigmaw,thl0bot,mu,mumin0,mumin1,mumin2,mulcl,mulclstar,winv,wlcl,ufrcinv,ufrclcl
    real(r8) :: thvj,exql,exqi,PGFc,rle,rkm,rpen,thl0top
    real(r8) :: qt0top,thl0lcl,qt0lcl,thv0lcl,thv0rel,rho0inv,thvubot,thvutop,autodet
    real(r8) :: rkfre,thv0j,rho0j,tj,rdqsdt
    real(r8) :: rbeta,ths0,aquad,bquad,cquad,xc1,xc2,excessu,excess0,xsat,xs,xs1,xs2
    real(r8) :: bogbot,bogtop,delbog,expfac,expfach,rhos0j,ppen,rcwp,rlwp,riwp,ufrcbelow,rainflx,ppen2
    real(r8) :: dpnb,delbognb,dragenb,expfacnb
    real(r8) :: snowflx,rmaxfrac,drage,qtdef,qpdef(ncnst),qcubelow,qlubelow,qiubelow,criqc
    real(r8) :: epsvarw                                                     ! vertical velocity variance at inversion base by meso-scale component.  
    real(r8) :: es(1)                                                       !  saturation vapor pressure
    real(r8) :: qs(1)                                                       !  saturation spec. humidity
    real(r8) :: gam(1)                                                      !  (L/cp)*dqs/dT
    integer  :: status                                                      !  return status of qsat call
    real(r8) :: thvlmin                                                     !  minimum theta_vl in the PBL

    real(r8)    erfarg,qsat_arg,thvlsrc             
    real(r8)    dpi                                                         ! Full(for internal) or half(for external) interface thickness for tkeavg
    real(r8)    x1,x2,f1,f2,xmid,f_thl,fmid,dx,j                            ! Used for 'thlsrc' calculation from ' p, qtsrc, thvmin' 
    integer     kthvmin,iter_scaleh,iter_xc  
    integer     jj
    real(r8) :: xsrc, xmean, xtop, xbot, xflx(0:mkx)
    real(r8) :: qt_dt, thl_dt

    !----- Some diagnostic internal output variables

    real(r8)  ufrcinvbase_out(mix)  !  Cumulus updraft fraction at PBL top
    real(r8)  ufrclcl_out(mix)      !  Cumulus updraft fraction at LCL ( or PBL top when LCL is below PBL top )
    real(r8)  winvbase_out(mix)     !  Cumulus updraft velocity at PBL top
    real(r8)  wlcl_out(mix)         !  Cumulus updraft velocity at LCL ( or PBL top when LCL is below PBL top )
    real(r8)  plcl_out(mix)         !  LCL pressure of source air
    real(r8)  pinv_out(mix)         !  PBL top pressure
    real(r8)  plfc_out(mix)         !  LFC pressure of source air
    real(r8)  pbup_out(mix)         !  Highest level of Positive Buoyancy
    real(r8)  ppen_out(mix)         !  Highest Level where Cu W = 0
    real(r8)  qtsrc_out(mix)        !  Sourse air qt
    real(r8)  thlsrc_out(mix)       !  Sourse air thl
    real(r8)  thvlsrc_out(mix)      !  Sourse air thvl
    real(r8)  emfkbup_out(mix)      !  Penetrative mass flux at 'kbup' interface
    real(r8)  cinlclh_out(mix)      !  convective inhibition upto LCL (CIN) ([J/kg])
    real(r8)  tkeavg_out(mix)       !  tkeavg over the PBL
    real(r8)  cbmflimit_out(mix)    !  Cbmf limiter
    real(r8)  zinv_out(mix)         !  PBL top height
    real(r8)  rcwp_out(mix)         !  Layer mean Cumulus LWP+IWP [ kg/m2 ] 
    real(r8)  rlwp_out(mix)         !  Layer mean Cumulus LWP [ kg/m2 ] 
    real(r8)  riwp_out(mix)         !  Layer mean Cumulus IWP [ kg/m2 ] 
    real(r8)  wu_out(mix,0:mkx)     !  Cumulus updraft vertical velocity ( defined from the release level, i.e. LCL or PBL top,
                                    !                                      to 'kpen-1' interface. Zero for the other interfaces) 
    real(r8)  qtu_out(mix,0:mkx)
    real(r8)  thlu_out(mix,0:mkx)
    real(r8)  thvu_out(mix,0:mkx)
    real(r8)  uu_out(mix,0:mkx)
    real(r8)  vu_out(mix,0:mkx)
    real(r8)  qtu_emf_out(mix,0:mkx) 
    real(r8)  thlu_emf_out(mix,0:mkx)  
    real(r8)  uu_emf_out(mix,0:mkx)   
    real(r8)  vu_emf_out(mix,0:mkx)   
    real(r8)  uemf_out(mix,0:mkx)   !  Net upward cumulus mass flux including penetrative entrainment flux (umf+emf)

    real(r8)  wu_s(0:mkx)           !  Cumulus updraft vertical velocity
    real(r8)  qtu_s(0:mkx)
    real(r8)  thlu_s(0:mkx)
    real(r8)  thvu_s(0:mkx)
    real(r8)  uu_s(0:mkx)
    real(r8)  vu_s(0:mkx)
    real(r8)  qtu_emf_s(0:mkx) 
    real(r8)  thlu_emf_s(0:mkx)  
    real(r8)  uu_emf_s(0:mkx)   
    real(r8)  vu_emf_s(0:mkx)
    real(r8)  uemf_s(0:mkx)   

    real(r8)  dwten_out(mix,mkx)
    real(r8)  diten_out(mix,mkx)
    real(r8)  flxrain_out(mix,0:mkx)  
    real(r8)  flxsnow_out(mix,0:mkx)  
    real(r8)  ntraprd_out(mix,mkx)    
    real(r8)  ntsnprd_out(mix,mkx)    

    real(r8)  dwten_s(mkx)
    real(r8)  diten_s(mkx)
    real(r8)  flxrain_s(0:mkx)  
    real(r8)  flxsnow_s(0:mkx)  
    real(r8)  ntraprd_s(mkx)    
    real(r8)  ntsnprd_s(mkx)    

    real(r8)  excessu_arr_out(mix,mkx)
    real(r8)  excessu_arr(mkx) 
    real(r8)  excessu_arr_s(mkx)
    real(r8)  excess0_arr_out(mix,mkx)
    real(r8)  excess0_arr(mkx)
    real(r8)  excess0_arr_s(mkx)
    real(r8)  xc_arr_out(mix,mkx)
    real(r8)  xc_arr(mkx)
    real(r8)  xc_arr_s(mkx)
    real(r8)  aquad_arr_out(mix,mkx)
    real(r8)  aquad_arr(mkx)
    real(r8)  aquad_arr_s(mkx)
    real(r8)  bquad_arr_out(mix,mkx)
    real(r8)  bquad_arr(mkx)
    real(r8)  bquad_arr_s(mkx)
    real(r8)  cquad_arr_out(mix,mkx) 
    real(r8)  cquad_arr(mkx)
    real(r8)  cquad_arr_s(mkx)
    real(r8)  bogbot_arr_out(mix,mkx)
    real(r8)  bogbot_arr(mkx)
    real(r8)  bogbot_arr_s(mkx)
    real(r8)  bogtop_arr_out(mix,mkx)
    real(r8)  bogtop_arr(mkx)
    real(r8)  bogtop_arr_s(mkx)

    real(r8)  exit_UWCu(mix)
    real(r8)  exit_conden(mix)
    real(r8)  exit_klclmkx(mix)
    real(r8)  exit_klfcmkx(mix)
    real(r8)  exit_ufrc(mix)
    real(r8)  exit_wtw(mix)
    real(r8)  exit_drycore(mix)
    real(r8)  exit_wu(mix)
    real(r8)  exit_cufilter(mix)
    real(r8)  exit_kinv1(mix)
    real(r8)  exit_rei(mix)

    real(r8)  limit_shcu(mix)
    real(r8)  limit_negcon(mix)
    real(r8)  limit_ufrc(mix)
    real(r8)  limit_ppen(mix)
    real(r8)  limit_emf(mix)
    real(r8)  limit_cinlcl(mix)
    real(r8)  limit_cin(mix)
    real(r8)  limit_cbmf(mix)
    real(r8)  limit_rei(mix)
    real(r8)  ind_delcin(mix)

    real(r8) :: ufrcinvbase_s, ufrclcl_s, winvbase_s, wlcl_s, plcl_s, pinv_s, plfc_s, &
                qtsrc_s, thlsrc_s, thvlsrc_s, emfkbup_s, cinlcl_s, pbup_s, ppen_s, cbmflimit_s, &
                tkeavg_s, zinv_s, rcwp_s, rlwp_s, riwp_s 
    real(r8) :: ufrcinvbase, winvbase, pinv, zinv, emfkbup, cbmflimit, rho0rel  

    !----- Variables for implicit CIN computation

    real(r8), dimension(mkx)         :: qv0_s  , ql0_s   , qi0_s   , s0_s    , u0_s    ,           & 
                                        v0_s   , t0_s    , qt0_s   , thl0_s  , thvl0_s , qvten_s , &
                                        qlten_s, qiten_s , qrten_s , qsten_s , sten_s  ,           &
                                        uten_s , vten_s  , cufrc_s , qcu_s   , qlu_s   , qiu_s   , &
                                        fer_s  , fdr_s   , qc_s    , qtten_s , slten_s 
    real(r8), dimension(0:mkx)       :: umf_s  , slflx_s , qtflx_s , ufrc_s  , uflx_s , vflx_s
    real(r8)                         :: cush_s , precip_s, snow_s  , cin_s   , rliq_s, cbmf_s, cnt_s, cnb_s
    real(r8)                         :: cin_i,cin_f,del_CIN,ke,alpha,thlj
    real(r8)                         :: cinlcl_i,cinlcl_f,del_cinlcl
    integer                          :: iter

    !----- Variables for temporary storages

    real(r8), dimension(mkx)   :: qv0_o, ql0_o, qi0_o, t0_o, s0_o, u0_o, v0_o
    real(r8), dimension(mkx)   :: qt0_o    , thl0_o   , thvl0_o   ,                         &
                                  qvten_o  , qlten_o  , qiten_o   , qrten_o   , qsten_o ,   &
                                  sten_o   , uten_o   , vten_o    , qcu_o     , qlu_o   ,   & 
                                  qiu_o    , cufrc_o  ,                                     &
                                  thv0bot_o, thv0top_o, thvl0bot_o, thvl0top_o,             &
                                  ssthl0_o , ssqt0_o  , ssu0_o    , ssv0_o    , qc_o    ,   &
                                  qtten_o  , slten_o  
    real(r8), dimension(0:mkx) :: umf_o    , slflx_o  , qtflx_o   , ufrc_o 
    real(r8), dimension(mix)   :: cush_o   , precip_o , snow_o    , rliq_o, cbmf_o, cnt_o, cnb_o
    real(r8), dimension(0:mkx) :: uflx_o   , vflx_o
    real(r8)                   :: tkeavg_o , thvlmin_o, qtsrc_o  , thvlsrc_o, thlsrc_o ,    &
                                  usrc_o   , vsrc_o   , plcl_o   , plfc_o   ,               &
                                  thv0lcl_o, cinlcl_o 
    integer                    :: kmin_o   , kinv_o   , klcl_o   , klfc_o  

    ! ------------------ !
    !                    !
    ! Define Parameters  !
    !                    !
    ! ------------------ !

    ! --------------------------------------------------------------------- !
    ! Choose new (correct saturated+unsaturated) or old (wrong unsaturated) !
    ! buoyancy sorting scheme.                                              !
    ! --------------------------------------------------------------------- !

    logical , parameter              :: use_newbuosort = .true.

    ! --------------------------------------------------------------------- !
    ! Set 'fer(kpen) = zero' (.true.) or use explicitly calculated non-zero !
    ! fer(kpen) (.false.)                                                   !
    ! --------------------------------------------------------------------- !

    logical , parameter              :: set_zeroferkpen = .false.

    ! ------------------------ !
    ! Iterative xc calculation !
    ! ------------------------ !

    integer , parameter              :: niter_xc = 3

    ! ----------------------------------------------------------- !
    ! Choice of 'CIN = cin' (.true.) or 'CIN = cinlcl' (.false.). !
    ! ----------------------------------------------------------- !

    logical , parameter              :: use_CINcin = .true.

    ! --------------------------------------------------------------- !
    ! Choice of 'explicit' ( 1 ) or 'implicit' ( 2 )  CIN.            !
    !                                                                 !
    ! When choose 'CIN = cinlcl' above,  it is recommended not to use ! 
    ! implicit CIN, i.e., do 'NOT' choose simultaneously :            !
    !            [ 'use_CINcin=.false. & 'iter_cin=2' ]               !
    ! since 'cinlcl' will be always set to zero whenever LCL is below !
    ! the PBL top interface in the current code. So, averaging cinlcl !
    ! of two iter_cin steps is likely not so good. Except that,   all !
    ! the other combinations of  'use_CINcin'  & 'iter_cin' are OK.   !
    !                                                                 !
    ! Feb 2007, Bundy: Note that use_CINcin = .false. will try to use !
    !           a variable (del_cinlcl) that is not currently set     !
    !                                                                 !
    ! --------------------------------------------------------------- !

    integer , parameter              :: iter_cin = 2

    ! ---------------------------------------------------------------- !
    ! Choice of 'self-detrainment' by negative buoyancy in calculating !
    ! cumulus updraft mass flux at the top interface in each layer.    !
    ! ---------------------------------------------------------------- !

    logical , parameter              :: use_self_detrain = .false.
    
    ! --------------------------------------------------------- !
    ! Cumulus momentum flux : turn-on (.true.) or off (.false.) !
    ! --------------------------------------------------------- !

    logical , parameter              :: use_momenflx = .true.

    ! ----------------------- !
    ! For lateral entrainment !
    ! ----------------------- !

    parameter (rle = 0.1_r8)         !  For critical stopping distance for lateral entrainment [no unit]
    parameter (rkm = 16.0_r8)        !  Determine the amount of air that is involved in buoyancy-sorting [no unit] 
    parameter (rpen = 10.0_r8)       !  For penetrative entrainment efficiency
    parameter (rkfre = 1.0_r8)       !  Vertical velocity variance as fraction of  tke. 
    parameter (rmaxfrac = 0.05_r8)   !  Maximum allowable updraft fraction
    parameter (mumin1 = 1.163_r8)    !  Normalized CIN ('mu') corresponding to 'rmaxfrac' at the PBL top
                                  !  obtaind by inverting 'rmaxfrac = 0.5*erfc(mumin1)'.
                                  !  [ rmaxfrac:mumin1 ] = [ 0.05:1.163, 0.1:0.906, 0.15:0.733, 0.2:0.595, 0.25:0.477 ] 
    parameter (rbuoy = 1.0_r8)       !  For nonhydrostatic pressure effects on updraft [no unit]
    parameter (rdrag = 1.0_r8)       !  Drag coefficient [no unit]

    parameter (epsvarw = 5.e-4_r8)   !  Variance of w at PBL top by meso-scale component [m2/s2]          
    parameter (PGFc = 0.7_r8)        !  This is used for calculating vertical variations cumulus  
                                  !  'u' & 'v' by horizontal PGF during upward motion [no unit]

    ! ---------------------------------------- !
    ! Bulk microphysics controlling parameters !
    ! --------------------------------------------------------------------------- ! 
    ! criqc    : Maximum condensate that can be hold by cumulus updraft [kg/kg]   !
    ! frc_rasn : Fraction of precipitable condensate in the expelled cloud water  !
    !            from cumulus updraft. The remaining fraction ('1-frc_rasn')  is  !
    !            'suspended condensate'.                                          !
    !                0 : all expelled condensate is 'suspended condensate'        ! 
    !                1 : all expelled condensate is 'precipitable condensate'     !
    ! kevp     : Evaporative efficiency                                           !
    ! noevap_krelkpen : No evaporation from 'krel' to 'kpen' layers               ! 
    ! --------------------------------------------------------------------------- !    

    parameter ( criqc = 1.e-3_r8 )        
    parameter ( frc_rasn = 1.0_r8 )
    parameter ( kevp = 2.e-6_r8 )
    logical , parameter :: noevap_krelkpen = .false.

    !------------------------!
    !                        !
    ! Start Main Calculation !
    !                        !
    !------------------------!

    ! ------------------------------------------------------- !
    ! Initialize output variables defined for all grid points !
    ! ------------------------------------------------------- !

    umf_out(:iend,0:mkx)         = 0.0_r8
    slflx_out(:iend,0:mkx)       = 0.0_r8
    qtflx_out(:iend,0:mkx)       = 0.0_r8
    qvten_out(:iend,:mkx)        = 0.0_r8
    qlten_out(:iend,:mkx)        = 0.0_r8
    qiten_out(:iend,:mkx)        = 0.0_r8
    sten_out(:iend,:mkx)         = 0.0_r8
    uten_out(:iend,:mkx)         = 0.0_r8
    vten_out(:iend,:mkx)         = 0.0_r8
    qrten_out(:iend,:mkx)        = 0.0_r8
    qsten_out(:iend,:mkx)        = 0.0_r8
    precip_out(:iend)            = 0.0_r8
    snow_out(:iend)              = 0.0_r8
    cufrc_out(:iend,:mkx)        = 0.0_r8
    qcu_out(:iend,:mkx)          = 0.0_r8
    qlu_out(:iend,:mkx)          = 0.0_r8
    qiu_out(:iend,:mkx)          = 0.0_r8
    fer_out(:iend,:mkx)          = 0.0_r8
    fdr_out(:iend,:mkx)          = 0.0_r8
    cinh_out(:iend)              = -1.0_r8
    cinlclh_out(:iend)           = -1.0_r8
    cbmf_out(:iend)              = 0.0_r8
    qc_out(:iend,:mkx)           = 0.0_r8
    rliq_out(:iend)              = 0.0_r8
    cnt_out(:iend)               = real(mkx, r8)
    cnb_out(:iend)               = 0.0_r8
    qtten_out(:iend,:mkx)        = 0.0_r8
    slten_out(:iend,:mkx)        = 0.0_r8
    ufrc_out(:iend,0:mkx)        = 0.0_r8

    uflx_out(:iend,0:mkx)        = 0.0_r8
    vflx_out(:iend,0:mkx)        = 0.0_r8
    
    ufrcinvbase_out(:iend)       = 0.0_r8
    ufrclcl_out(:iend)           = 0.0_r8
    winvbase_out(:iend)          = 0.0_r8
    wlcl_out(:iend)              = 0.0_r8
    plcl_out(:iend)              = 0.0_r8
    pinv_out(:iend)              = 0.0_r8
    plfc_out(:iend)              = 0.0_r8
    pbup_out(:iend)              = 0.0_r8
    ppen_out(:iend)              = 0.0_r8
    qtsrc_out(:iend)             = 0.0_r8
    thlsrc_out(:iend)            = 0.0_r8
    thvlsrc_out(:iend)           = 0.0_r8
    emfkbup_out(:iend)           = 0.0_r8
    cbmflimit_out(:iend)         = 0.0_r8
    tkeavg_out(:iend)            = 0.0_r8
    zinv_out(:iend)              = 0.0_r8
    rcwp_out(:iend)              = 0.0_r8
    rlwp_out(:iend)              = 0.0_r8
    riwp_out(:iend)              = 0.0_r8

    wu_out(:iend,0:mkx)          = 0.0_r8
    qtu_out(:iend,0:mkx)         = 0.0_r8
    thlu_out(:iend,0:mkx)        = 0.0_r8
    thvu_out(:iend,0:mkx)        = 0.0_r8
    uu_out(:iend,0:mkx)          = 0.0_r8
    vu_out(:iend,0:mkx)          = 0.0_r8
    qtu_emf_out(:iend,0:mkx)     = 0.0_r8
    thlu_emf_out(:iend,0:mkx)    = 0.0_r8
    uu_emf_out(:iend,0:mkx)      = 0.0_r8
    vu_emf_out(:iend,0:mkx)      = 0.0_r8
    uemf_out(:iend,0:mkx)        = 0.0_r8

    dwten_out(:iend,:mkx)        = 0.0_r8
    diten_out(:iend,:mkx)        = 0.0_r8
    flxrain_out(:iend,0:mkx)     = 0.0_r8  
    flxsnow_out(:iend,0:mkx)     = 0.0_r8
    ntraprd_out(:iend,mkx)       = 0.0_r8
    ntsnprd_out(:iend,mkx)       = 0.0_r8

    excessu_arr_out(:iend,:mkx)  = 0.0_r8
    excess0_arr_out(:iend,:mkx)  = 0.0_r8
    xc_arr_out(:iend,:mkx)       = 0.0_r8
    aquad_arr_out(:iend,:mkx)    = 0.0_r8
    bquad_arr_out(:iend,:mkx)    = 0.0_r8
    cquad_arr_out(:iend,:mkx)    = 0.0_r8
    bogbot_arr_out(:iend,:mkx)   = 0.0_r8
    bogtop_arr_out(:iend,:mkx)   = 0.0_r8

    exit_UWCu(:iend)             = 0.0_r8 
    exit_conden(:iend)           = 0.0_r8 
    exit_klclmkx(:iend)          = 0.0_r8 
    exit_klfcmkx(:iend)          = 0.0_r8 
    exit_ufrc(:iend)             = 0.0_r8 
    exit_wtw(:iend)              = 0.0_r8 
    exit_drycore(:iend)          = 0.0_r8 
    exit_wu(:iend)               = 0.0_r8 
    exit_cufilter(:iend)         = 0.0_r8 
    exit_kinv1(:iend)            = 0.0_r8 
    exit_rei(:iend)              = 0.0_r8 

    limit_shcu(:iend)            = 0.0_r8 
    limit_negcon(:iend)          = 0.0_r8 
    limit_ufrc(:iend)            = 0.0_r8
    limit_ppen(:iend)            = 0.0_r8
    limit_emf(:iend)             = 0.0_r8
    limit_cinlcl(:iend)          = 0.0_r8
    limit_cin(:iend)             = 0.0_r8
    limit_cbmf(:iend)            = 0.0_r8
    limit_rei(:iend)             = 0.0_r8

    ind_delcin(:iend)            = 0.0_r8

    !---------------------------------------------------------!
    !                                                         !
    ! Start the big i loop where i is a horozontal grid index !
    !                                                         !
    !---------------------------------------------------------!

    do i = 1, iend                                        ! start of big i loop

      id_exit = .false.

      ! --------------------------------------------- !
      ! Define new input variables at each grid point !
      ! --------------------------------------------- !

      ps0(0:mkx)       = ps0_in(i,0:mkx)
      zs0(0:mkx)       = zs0_in(i,0:mkx)
      p0(:mkx)         = p0_in(i,:mkx)
      z0(:mkx)         = z0_in(i,:mkx)
      dp0(:mkx)        = dp0_in(i,:mkx)
      u0(:mkx)         = u0_in(i,:mkx)
      v0(:mkx)         = v0_in(i,:mkx)
      qv0(:mkx)        = qv0_in(i,:mkx)
      ql0(:mkx)        = ql0_in(i,:mkx)
      qi0(:mkx)        = qi0_in(i,:mkx)
      t0(:mkx)         = t0_in(i,:mkx)
      s0(:mkx)         = s0_in(i,:mkx)
      tke(0:mkx)       = tke_in(i,0:mkx)
      cldfrct(:mkx)    = cldfrct_in(i,:mkx)
      concldfrct(:mkx) = concldfrct_in(i,:mkx)
      pblh             = pblh_in(i)
      cush             = cush_inout(i)

      ! --------------------------------------------------------- !
      ! Compute other basic thermodynamic variables directly from ! 
      ! the input variables at each grid point                    !
      ! --------------------------------------------------------- !

      !----- 1. Compute internal environmental variables
      
      exn0(:mkx)   = (p0(:mkx)/p00)**rovcp
      exns0(0:mkx) = (ps0(0:mkx)/p00)**rovcp
      qt0(:mkx)    = (qv0(:mkx) + ql0(:mkx) + qi0(:mkx))
      thl0(:mkx)   = (t0(:mkx) - xlv*ql0(:mkx)/cp - xls*qi0(:mkx)/cp)/exn0(:mkx)
      thvl0(:mkx)  = (1._r8 + zvir*qt0(:mkx))*thl0(:mkx)

      !----- 2. Compute slopes of environmental variables

      ssthl0       = slope(mkx,thl0,p0) ! Dimension of ssthl0(:mkx) is implicit
      ssqt0        = slope(mkx,qt0 ,p0)
      ssu0         = slope(mkx,u0  ,p0)
      ssv0         = slope(mkx,v0  ,p0)

      !----- 3. Compute "thv0" and "thvl0" at top/bottom interfaces

      do k = 1, mkx

        thl0bot = thl0(k) + ssthl0(k)*(ps0(k-1) - p0(k))
        qt0bot  = qt0(k)  + ssqt0(k) *(ps0(k-1) - p0(k))
        call conden(ps0(k-1),thl0bot,qt0bot,thj,qvj,qlj,qij,qse,id_check,qsat)
        if(id_check.eq.1) then
          exit_conden(i) = 1._r8
          id_exit = .true.
          go to 333
        end if
        thv0bot(k)  = thj*(1._r8 + zvir*qvj - qlj - qij)
        thvl0bot(k) = thl0bot*(1._r8 + zvir*qt0bot)
          
        thl0top       = thl0(k) + ssthl0(k)*(ps0(k) - p0(k))
        qt0top        =  qt0(k) + ssqt0(k) *(ps0(k) - p0(k))
        call conden(ps0(k),thl0top,qt0top,thj,qvj,qlj,qij,qse,id_check,qsat)
        if(id_check.eq.1) then
          exit_conden(i) = 1._r8
          id_exit = .true.
          go to 333
        end if 
        thv0top(k)  = thj*(1._r8 + zvir*qvj - qlj - qij)
        thvl0top(k) = thl0top*(1._r8 + zvir*qt0top)

      end do

      !----- 4. Compute "qs0(k)" at layer mid-point for condensate sink calculation

      do k = 1, mkx
        call conden(p0(k),thl0(k),qt0(k),thj,qvj,qlj,qij,qse,id_check,qsat)
        if(id_check.eq.1) then
          exit_conden(i) = 1._r8
          id_exit = .true.
          go to 333
        end if 
        qs0(k) = qse
      end do 

      ! -------------------------------------------------------------------- !
      ! Save input and related thermodynamic variables ( not associated with !
      ! cumulus convection ) for use at "iter_cin=2" when "del_CIN >= 0"     !
      ! -------------------------------------------------------------------- !

      qv0_o(:mkx)          = qv0(:mkx)
      ql0_o(:mkx)          = ql0(:mkx)
      qi0_o(:mkx)          = qi0(:mkx)
      t0_o(:mkx)           = t0(:mkx)
      s0_o(:mkx)           = s0(:mkx)
      u0_o(:mkx)           = u0(:mkx)
      v0_o(:mkx)           = v0(:mkx)
      qt0_o(:mkx)          = qt0(:mkx)
      thl0_o(:mkx)         = thl0(:mkx)
      thvl0_o(:mkx)        = thvl0(:mkx)
      ssthl0_o(:mkx)       = ssthl0(:mkx)
      ssqt0_o(:mkx)        = ssqt0(:mkx)
      thv0bot_o(:mkx)      = thv0bot(:mkx)
      thv0top_o(:mkx)      = thv0top(:mkx)
      thvl0bot_o(:mkx)     = thvl0bot(:mkx)
      thvl0top_o(:mkx)     = thvl0top(:mkx)
      ssu0_o(:mkx)         = ssu0(:mkx) 
      ssv0_o(:mkx)         = ssv0(:mkx) 

      ! ------------------------------------------------------------------------ !
      ! Initialize output variables before cumulus convection at each grid point !
      ! ------------------------------------------------------------------------ !

      umf(0:mkx)          = 0.0_r8
      emf(0:mkx)          = 0.0_r8
      slflx(0:mkx)        = 0.0_r8
      qtflx(0:mkx)        = 0.0_r8
      uflx(0:mkx)         = 0.0_r8
      vflx(0:mkx)         = 0.0_r8
      qvten(:mkx)         = 0.0_r8
      qlten(:mkx)         = 0.0_r8
      qiten(:mkx)         = 0.0_r8
      sten(:mkx)          = 0.0_r8
      uten(:mkx)          = 0.0_r8
      vten(:mkx)          = 0.0_r8
      qrten(:mkx)         = 0.0_r8
      qsten(:mkx)         = 0.0_r8
      dwten(:mkx)         = 0.0_r8
      diten(:mkx)         = 0.0_r8
      precip              = 0.0_r8
      snow                = 0.0_r8
      cufrc(:mkx)         = 0.0_r8
      qcu(:mkx)           = 0.0_r8
      qlu(:mkx)           = 0.0_r8
      qiu(:mkx)           = 0.0_r8
      fer(:mkx)           = 0.0_r8
      fdr(:mkx)           = 0.0_r8
      cin                 = 0.0_r8
      cbmf                = 0.0_r8
      qc(:mkx)            = 0.0_r8
      qc_l(:mkx)          = 0.0_r8
      qc_i(:mkx)          = 0.0_r8
      rliq                = 0.0_r8
      cnt                 = real(mkx, r8)
      cnb                 = 0.0_r8
      qtten(:mkx)         = 0.0_r8
      slten(:mkx)         = 0.0_r8   
      ufrc(0:mkx)         = 0.0_r8  

      thlu(0:mkx)         = 0.0_r8
      qtu(0:mkx)          = 0.0_r8
      uu(0:mkx)           = 0.0_r8
      vu(0:mkx)           = 0.0_r8
      wu(0:mkx)           = 0.0_r8
      thvu(0:mkx)         = 0.0_r8
      thlu_emf(0:mkx)     = 0.0_r8
      qtu_emf(0:mkx)      = 0.0_r8
      uu_emf(0:mkx)       = 0.0_r8
      vu_emf(0:mkx)       = 0.0_r8
      
      ufrcinvbase         = 0.0_r8
      ufrclcl             = 0.0_r8
      winvbase            = 0.0_r8
      wlcl                = 0.0_r8
      emfkbup             = 0.0_r8 
      cbmflimit           = 0.0_r8
      excessu_arr(:mkx)   = 0.0_r8
      excess0_arr(:mkx)   = 0.0_r8
      xc_arr(:mkx)        = 0.0_r8
      aquad_arr(:mkx)     = 0.0_r8
      bquad_arr(:mkx)     = 0.0_r8
      cquad_arr(:mkx)     = 0.0_r8
      bogbot_arr(:mkx)    = 0.0_r8
      bogtop_arr(:mkx)    = 0.0_r8

      uemf(0:mkx)         = 0.0_r8
      comsub(:mkx)        = 0.0_r8
      qcten_sink(:mkx)    = 0.0_r8
      qlten_sink(:mkx)    = 0.0_r8
      qiten_sink(:mkx)    = 0.0_r8 

    !-----------------------------------------------! 
    ! Below 'iter' loop is for implicit CIN closure !
    !-----------------------------------------------!

    ! ----------------------------------------------------------------------------- ! 
    ! It is important to note that this iterative cin loop is located at the outest !
    ! shell of the code. Thus, source air properties can also be changed during the !
    ! iterative cin calculation, because cumulus convection induces non-zero fluxes !
    ! even at interfaces below PBL top height through 'fluxbelowinv' subroutine.    !
    ! ----------------------------------------------------------------------------- !

    do iter = 1, iter_cin

       ! ---------------------------------------------------------------------- ! 
       ! Cumulus scale height                                                   ! 
       ! In contrast to the premitive code, cumulus scale height is iteratively !
       ! calculated at each time step, and at each iterative cin step.          !
       ! It is not clear whether I should locate below two lines within or  out !
       ! of the iterative cin loop.                                             !
       ! ---------------------------------------------------------------------- !

       tscaleh = cush                        
       cush    = -1._r8

       ! ----------------------------------------------------------------------- !
       ! Find PBL top height interface index, 'kinv-1' where 'kinv' is the layer !
       ! index with PBLH in it. When PBLH is exactly at interface, 'kinv' is the !
       ! layer index having PBLH as a lower interface.                           !
       ! In the previous code, I set the lower limit of 'kinv' by 2  in order to !
       ! be consistent with the other parts of the code. However in the modified !
       ! code, I allowed 'kinv' to be 1 & if 'kinv = 1', I just exit the program !
       ! without performing cumulus convection. This new approach seems to be    !
       ! more reasonable: if PBL height is within 'kinv=1' layer, surface is STL !
       ! interface (bflxs <= 0) and interface just above the surface should be   !
       ! either non-turbulent (Ri>0.19) or stably turbulent (0<=Ri<0.19 but this !
       ! interface is identified as a base external interface of upperlying CL.  !
       ! Thus, when 'kinv=1', PBL scheme guarantees 'bflxs <= 0'.  For this case !
       ! it is reasonable to assume that cumulus convection does not happen.     !
       ! When these is SBCL, PBL height from the PBL scheme is likely to be very !
       ! close at 'kinv-1' interface, but not exactly, since 'zi' information is !
       ! changed between two model time steps. In order to ensure correct identi !
       ! fication of 'kinv' for general case including SBCL, I imposed an offset !
       ! of 5 [m] in the below 'kinv' finding block.                             !
       ! ----------------------------------------------------------------------- !
       
       do k = mkx - 1, 1, -1 
         if((pblh + 5._r8 - zs0(k))*(pblh + 5._r8 - zs0(k+1)) .lt. 0._r8) then
           kinv = k + 1 
           go to 15
         endif 
       end do
       kinv = 1
15     continue    

       if(kinv.le.1) then          
        ! write(iulog,*) 'PBL height is within the lowest model layer'
          exit_kinv1(i) = 1._r8
          id_exit = .true.
          go to 333
       endif
       ! From here, it must be 'kinv >= 2' from now on. 

       ! -------------------------------------------------------------------------- !
       ! Find PBL averaged tke ('tkeavg') and minimum 'thvl' ('thvlmin') in the PBL !
       ! In the current code, 'tkeavg' is obtained by averaging all interfacial TKE !
       ! within the PBL. However, in order to be conceptually consistent with   PBL !
       ! scheme, 'tkeavg' should be calculated by considering surface buoyancy flux.!
       ! If surface buoyancy flux is positive ( bflxs >0 ), surface interfacial TKE !
       ! should be included in calculating 'tkeavg', while if bflxs <= 0,   surface !
       ! interfacial TKE should not be included in calculating 'tkeavg'.   I should !
       ! modify the code when 'bflxs' is available as an input of cumulus scheme.   !
       ! 'thvlmin' is a minimum 'thvl' within PBL obtained by comparing top &  base !
       ! interface values of 'thvl' in each layers within the PBL.                  !
       ! -------------------------------------------------------------------------- !
       
       dpsum    = 0._r8
       tkeavg   = 0._r8
       thvlmin  = 1000._r8
       kmin     = 1
       do k = 0, kinv - 1   ! Here, 'k' is an interfacial layer index.  
         if(k .eq. 0) then
            dpi = ps0(0) - p0(1)
         elseif(k .eq. (kinv-1)) then 
            dpi = p0(kinv-1) - ps0(kinv-1)
         else
            dpi = p0(k) - p0(k+1)
         endif 
         dpsum  = dpsum  + dpi  
         tkeavg = tkeavg + dpi*tke(k) 
         if( k.ne.0 ) thvlmin = min(thvlmin,min(thvl0bot(k),thvl0top(k)))
       end do
       tkeavg  = tkeavg/dpsum

       ! ------------------------------------------------------------------ !
       ! Find characteristics of cumulus source air: qtsrc,thlsrc,usrc,vsrc !
       ! Note that 'thlsrc' was con-cocked using 'thvlsrc' and 'qtsrc'.     !
       ! 'qtsrc' is defined as the lowest layer mid-point value;   'thlsrc' !
       ! is from 'qtsrc' and 'thvlmin=thvlsrc'; 'usrc' & 'vsrc' are defined !
       ! as the values just below the PBL top interface.                    !
       ! ------------------------------------------------------------------ !

       qtsrc   = qt0(1)                     
       thvlsrc = thvlmin 
       thlsrc  = thvlsrc / ( 1._r8 + zvir * qtsrc )  
       usrc    = u0(kinv-1) + ssu0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )             
       vsrc    = v0(kinv-1) + ssv0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )             

       ! ------------------------------------------------------------------ !
       ! Find LCL of the source air and a layer index containing LCL (klcl) !
       ! When the LCL is exactly at the interface, 'klcl' is a layer index  ! 
       ! having 'plcl' as the lower interface similar to the 'kinv' case.   !
       ! In the previous code, I assumed that if LCL is located within the  !
       ! lowest model layer ( 1 ) or the top model layer ( mkx ), then  no  !
       ! convective adjustment is performed and just exited.   However, in  !
       ! the revised code, I relaxed the first constraint and  even though  !
       ! LCL is at the lowest model layer, I allowed cumulus convection to  !
       ! be initiated. For this case, cumulus convection should be started  !
       ! from the PBL top height, as shown in the following code.           !
       ! When source air is already saturated even at the surface, klcl is  !
       ! set to 1.                                                          !
       ! ------------------------------------------------------------------ !

       plcl = qsinvert(qtsrc,thlsrc,ps0(0),qsat)
       do k = 0, mkx
         if( ps0(k) .lt. plcl ) then
           klcl = k
           go to 25
         endif           
       end do
       klcl = mkx
25     continue
       klcl = max(1,klcl)

       if(plcl.lt.30000.) then               
     ! if(klcl.eq.mkx) then          
        ! write(iulog,*) 'Too dry source air'
          exit_klclmkx(i) = 1._r8
          id_exit = .true.
          go to 333
       endif

       ! ------------------------------------------------------------- !
       ! Calculate environmental virtual potential temperature at LCL, !
       !'thv0lcl' which is solely used in the 'cin' calculation. Note  !
       ! that 'thv0lcl' is calculated first by calculating  'thl0lcl'  !
       ! and 'qt0lcl' at the LCL, and performing 'conden' afterward,   !
       ! in fully consistent with the other parts of the code.         !
       ! ------------------------------------------------------------- !

       thl0lcl = thl0(klcl) + ssthl0(klcl) * ( plcl - p0(klcl) )
       qt0lcl  = qt0(klcl)  + ssqt0(klcl)  * ( plcl - p0(klcl) )
       call conden(plcl,thl0lcl,qt0lcl,thj,qvj,qlj,qij,qse,id_check,qsat)
       if(id_check.eq.1) then
          exit_conden(i) = 1._r8
          id_exit = .true.
          go to 333
       end if
       thv0lcl = thj * (1._r8 + zvir * qvj - qlj - qij )

       ! ------------------------------------------------------------------------ !
       ! Compute Convective Inhibition, 'cin' & 'cinlcl' [J/kg]=[m2/s2] TKE unit. !
       !                                                                          !
       ! 'cin' (cinlcl) is computed from the PBL top interface to LFC (LCL) using ! 
       ! piecewisely reconstructed environmental profiles, assuming environmental !
       ! buoyancy profile within each layer ( or from LCL to upper interface in   !
       ! each layer ) is simply a linear profile. For the purpose of cin (cinlcl) !
       ! calculation, we simply assume that lateral entrainment does not occur in !
       ! updrafting cumulus plume, i.e., cumulus source air property is conserved.!
       ! Below explains some rules used in the calculations of cin (cinlcl).   In !
       ! general, both 'cin' and 'cinlcl' are calculated from a PBL top interface !
       ! to LCL and LFC, respectively :                                           !
       ! 1. If LCL is lower than the PBL height, cinlcl = 0 and cin is calculated !
       !    from PBL height to LFC.                                               !
       ! 2. If LCL is higher than PBL height,   'cinlcl' is calculated by summing !
       !    both positive and negative cloud buoyancy up to LCL using 'single_cin'!
       !    From the LCL to LFC, however, only negative cloud buoyancy is counted !
       !    to calculate final 'cin' upto LFC.                                    !
       ! 3. If either 'cin' or 'cinlcl' is negative, they are set to be zero.     !
       ! In the below code, 'klfc' is the layer index containing 'LFC' similar to !
       ! 'kinv' and 'klcl'.                                                       !
       ! ------------------------------------------------------------------------ !

        cin    = 0._r8
        cinlcl = 0._r8
        plfc   = 0._r8
        klfc   = mkx

        ! ------------------------------------------------------------------------- !
        ! Case 1. LCL height is higher than PBL interface ( 'pLCL <= ps0(kinv-1)' ) !
        ! ------------------------------------------------------------------------- !

        if( klcl .ge. kinv ) then

           do k = kinv, mkx - 1
             if( k .lt. klcl ) then
                thvubot = thvlsrc
                thvutop = thvlsrc  
                cin     = cin + single_cin(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop)
             elseif( k .eq. klcl ) then
                !----- Bottom to LCL
                thvubot = thvlsrc
                thvutop = thvlsrc
                cin     = cin + single_cin(ps0(k-1),thv0bot(k),plcl,thv0lcl,thvubot,thvutop)
                if(cin.lt.0) limit_cinlcl(i) = 1._r8
                cinlcl  = max(cin,0._r8)
                cin     = cinlcl
                !----- LCL to Top
                thvubot = thvlsrc
                call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
                if(id_check.eq.1) then
                   exit_conden(i) = 1._r8
                   id_exit = .true.
                   go to 333
                end if
                thvutop = thj*(1._r8 + zvir*qvj - qlj - qij )
                call getbuoy(plcl,thv0lcl,ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
                if(plfc .gt. 0._r8) then 
                   klfc = k 
                   go to 35
                end if
             else
                thvubot = thvutop
                call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
                if(id_check.eq.1) then
                   exit_conden(i) = 1._r8
                   id_exit = .true.
                   go to 333
                end if
                thvutop = thj*(1._r8 + zvir*qvj - qlj - qij )
                call getbuoy(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
                if(plfc .gt. 0._r8) then 
                   klfc = k
                   go to 35
                end if 
            endif
          end do        

       ! ----------------------------------------------------------------------- !
       ! Case 2. LCL height is lower than PBL interface ( 'pLCL > ps0(kinv-1)' ) !
       ! ----------------------------------------------------------------------- !

       else
          cinlcl = 0._r8 
          do k = kinv, mkx - 1
             call conden(ps0(k-1),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
             if(id_check.eq.1) then
                exit_conden(i) = 1._r8
                id_exit = .true.
                go to 333
             end if
             thvubot = thj*(1._r8 + zvir*qvj - qlj - qij )
             call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
             if(id_check.eq.1) then
                exit_conden(i) = 1._r8
                id_exit = .true.
                go to 333
             end if
             thvutop = thj*(1._r8 + zvir*qvj - qlj - qij )
             call getbuoy(ps0(k-1),thv0bot(k),ps0(k),thv0top(k),thvubot,thvutop,plfc,cin)
             if(plfc .gt. 0._r8) then 
                klfc = k
                go to 35
             end if 
         end do
       endif  ! End of CIN case selection

 35    continue
       if(cin.lt.0) limit_cin(i) = 1._r8
       cin = max(0._r8,cin)
       if( klfc .ge. mkx ) then
           klfc = mkx
         ! write(iulog,*) 'klfc >= mkx'
           exit_klfcmkx(i) = 1._r8
           id_exit = .true.
           go to 333
       endif

       ! ---------------------------------------------------------------------- !
       ! In order to calculate implicit 'cin' (or 'cinlcl'), save the initially !
       ! calculated 'cin' and 'cinlcl', and other related variables. These will !
       ! be restored after calculating implicit CIN.                            !
       ! ---------------------------------------------------------------------- !

       if(iter .eq. 1) then 
          cin_i       = cin
          cinlcl_i    = cinlcl
          ke          = rbuoy / ( rkfre * tkeavg + epsvarw ) 
          kmin_o      = kmin     
          kinv_o      = kinv     
          klcl_o      = klcl     
          klfc_o      = klfc    
          plcl_o      = plcl    
          plfc_o      = plfc     
          tkeavg_o    = tkeavg   
          thvlmin_o   = thvlmin
          qtsrc_o     = qtsrc    
          thvlsrc_o   = thvlsrc  
          thlsrc_o    = thlsrc
          usrc_o      = usrc     
          vsrc_o      = vsrc     
          thv0lcl_o   = thv0lcl  
       endif   

       ! -------------------------------------------------------------------------- !
       ! Calculate implicit 'cin' by averaging initial and final cins.    Note that !
       ! implicit CIN is adopted only when cumulus convection stabilized the system,!
       ! i.e., only when 'del_CIN >0'. If 'del_CIN<=0', just use explicit CIN. Note !
       ! also that since 'cinlcl' is set to zero whenever LCL is below the PBL top, !
       ! (see above CIN calculation part), the use of 'implicit CIN=cinlcl'  is not !
       ! good. Thus, when using implicit CIN, always try to only use 'implicit CIN= !
       ! cin', not 'implicit CIN=cinlcl'. However, both 'CIN=cin' and 'CIN=cinlcl'  !
       ! are good when using explicit CIN.                                          !
       ! -------------------------------------------------------------------------- !

       if(iter .ne. 1) then

          cin_f = cin
          cinlcl_f = cinlcl
          if(use_CINcin) then
             del_CIN = cin_f - cin_i
          else
             del_CIN = cinlcl_f - cinlcl_i
          endif

          if(del_CIN .gt. 0._r8) then

             ! -------------------------------------------------------------- ! 
             ! Calculate implicit 'cin' and 'cinlcl'. Note that when we chose !
             ! to use 'implicit CIN = cin', choose 'cinlcl = cinlcl_i' below: !
             ! because iterative CIN only aims to obtain implicit CIN,  once  !
             ! we obtained 'implicit CIN=cin', it is good to use the original !
             ! profiles information for all the other variables after that.   !
             ! Note 'cinlcl' will be explicitly used in calculating  'wlcl' & !
             ! 'ufrclcl' after calculating 'winv' & 'ufrcinv'  at the PBL top !
             ! interface later, after calculating 'cbmf'.                     !
             ! -------------------------------------------------------------- !
           
             alpha = compute_alpha( del_CIN, ke ) 
             cin   = cin_i + alpha * del_CIN
             if(use_CINcin) then
                cinlcl = cinlcl_i                 
             else
                !bundy: warning, del_cinlcl needs to be set if use_CINcin = .false.
                cinlcl = cinlcl_i + alpha * del_cinlcl   
             endif

             ! ----------------------------------------------------------------- !
             ! Restore the original values from the previous 'iter_cin' step (1) !
             ! to compute correct tendencies for (n+1) time step by implicit CIN !
             ! ----------------------------------------------------------------- !

             kmin      = kmin_o    
             kinv      = kinv_o     
             klcl      = klcl_o     
             klfc      = klfc_o    
             plcl      = plcl_o    
             plfc      = plfc_o     
             tkeavg    = tkeavg_o   
             thvlmin   = thvlmin_o
             qtsrc     = qtsrc_o    
             thvlsrc   = thvlsrc_o  
             thlsrc    = thlsrc_o
             usrc      = usrc_o     
             vsrc      = vsrc_o     
             thv0lcl   = thv0lcl_o  

             qv0(:mkx)            = qv0_o(:mkx)
             ql0(:mkx)            = ql0_o(:mkx)
             qi0(:mkx)            = qi0_o(:mkx)
             t0(:mkx)             = t0_o(:mkx)
             s0(:mkx)             = s0_o(:mkx)
             u0(:mkx)             = u0_o(:mkx)
             v0(:mkx)             = v0_o(:mkx)
             qt0(:mkx)            = qt0_o(:mkx)
             thl0(:mkx)           = thl0_o(:mkx)
             thvl0(:mkx)          = thvl0_o(:mkx)
             ssthl0(:mkx)         = ssthl0_o(:mkx)
             ssqt0(:mkx)          = ssqt0_o(:mkx)
             thv0bot(:mkx)        = thv0bot_o(:mkx)
             thv0top(:mkx)        = thv0top_o(:mkx)
             thvl0bot(:mkx)       = thvl0bot_o(:mkx)
             thvl0top(:mkx)       = thvl0top_o(:mkx)
             ssu0(:mkx)           = ssu0_o(:mkx) 
             ssv0(:mkx)           = ssv0_o(:mkx) 

             ! ------------------------------------------------------ !
             ! Initialize all fluxes, tendencies, and other variables ! 
             ! in association with cumulus convection.                !
             ! ------------------------------------------------------ ! 

             umf(0:mkx)          = 0.0_r8
             emf(0:mkx)          = 0.0_r8
             slflx(0:mkx)        = 0.0_r8
             qtflx(0:mkx)        = 0.0_r8
             uflx(0:mkx)         = 0.0_r8
             vflx(0:mkx)         = 0.0_r8
             qvten(:mkx)         = 0.0_r8
             qlten(:mkx)         = 0.0_r8
             qiten(:mkx)         = 0.0_r8
             sten(:mkx)          = 0.0_r8
             uten(:mkx)          = 0.0_r8
             vten(:mkx)          = 0.0_r8
             qrten(:mkx)         = 0.0_r8
             qsten(:mkx)         = 0.0_r8
             dwten(:mkx)         = 0.0_r8
             diten(:mkx)         = 0.0_r8
             precip              = 0.0_r8
             snow                = 0.0_r8
             cufrc(:mkx)         = 0.0_r8
             qcu(:mkx)           = 0.0_r8
             qlu(:mkx)           = 0.0_r8
             qiu(:mkx)           = 0.0_r8
             fer(:mkx)           = 0.0_r8
             fdr(:mkx)           = 0.0_r8
             qc(:mkx)            = 0.0_r8
             qc_l(:mkx)          = 0.0_r8
             qc_i(:mkx)          = 0.0_r8
             rliq                = 0.0_r8
             cbmf                = 0.0_r8
             cnt                 = real(mkx, r8)
             cnb                 = 0.0_r8
             qtten(:mkx)         = 0.0_r8
             slten(:mkx)         = 0.0_r8
             ufrc(0:mkx)         = 0.0_r8

             thlu(0:mkx)         = 0.0_r8
             qtu(0:mkx)          = 0.0_r8
             uu(0:mkx)           = 0.0_r8
             vu(0:mkx)           = 0.0_r8
             wu(0:mkx)           = 0.0_r8
             thvu(0:mkx)         = 0.0_r8
             thlu_emf(0:mkx)     = 0.0_r8
             qtu_emf(0:mkx)      = 0.0_r8
             uu_emf(0:mkx)       = 0.0_r8
             vu_emf(0:mkx)       = 0.0_r8             

             ! -------------------------------------------------- !
             ! Below are diagnostic output variables for detailed !
             ! analysis of cumulus scheme.                        !
             ! -------------------------------------------------- ! 

             ufrcinvbase         = 0.0_r8
             ufrclcl             = 0.0_r8
             winvbase            = 0.0_r8
             wlcl                = 0.0_r8
             emfkbup             = 0.0_r8 
             cbmflimit           = 0.0_r8
             excessu_arr(:mkx)   = 0.0_r8
             excess0_arr(:mkx)   = 0.0_r8
             xc_arr(:mkx)        = 0.0_r8
             aquad_arr(:mkx)     = 0.0_r8
             bquad_arr(:mkx)     = 0.0_r8
             cquad_arr(:mkx)     = 0.0_r8
             bogbot_arr(:mkx)    = 0.0_r8
             bogtop_arr(:mkx)    = 0.0_r8

          else ! When 'del_CIN < 0', use explicit CIN instead of implicit CIN.
           
             ! ----------------------------------------------------------- ! 
             ! Identifier showing whether explicit or implicit CIN is used !
             ! ----------------------------------------------------------- ! 

             ind_delcin(i) = 1._r8             
   
             ! --------------------------------------------------------- !
             ! Restore original output values of "iter_cin = 1" and exit !
             ! --------------------------------------------------------- !

             umf_out(i,0:mkx)         = umf_s(0:mkx)
             qvten_out(i,:mkx)        = qvten_s(:mkx)
             qlten_out(i,:mkx)        = qlten_s(:mkx)  
             qiten_out(i,:mkx)        = qiten_s(:mkx)
             sten_out(i,:mkx)         = sten_s(:mkx)
             uten_out(i,:mkx)         = uten_s(:mkx)  
             vten_out(i,:mkx)         = vten_s(:mkx)
             qrten_out(i,:mkx)        = qrten_s(:mkx)
             qsten_out(i,:mkx)        = qsten_s(:mkx)  
             precip_out(i)            = precip_s
             snow_out(i)              = snow_s
             cush_inout(i)            = cush_s
             cufrc_out(i,:mkx)        = cufrc_s(:mkx)  
             slflx_out(i,0:mkx)       = slflx_s(0:mkx)  
             qtflx_out(i,0:mkx)       = qtflx_s(0:mkx)
             qcu_out(i,:mkx)          = qcu_s(:mkx)    
             qlu_out(i,:mkx)          = qlu_s(:mkx)  
             qiu_out(i,:mkx)          = qiu_s(:mkx)  
             cbmf_out(i)              = cbmf_s
             qc_out(i,:mkx)           = qc_s(:mkx)  
             rliq_out(i)              = rliq_s
             cnt_out(i)               = cnt_s
             cnb_out(i)               = cnb_s
             
             ! -------------------------------------------------- ! 
             ! Below are diagnostic output variables for detailed !
             ! analysis of cumulus scheme.                        !
             ! -------------------------------------------------- !   

             fer_out(i,mkx:1:-1)      = fer_s(:mkx)  
             fdr_out(i,mkx:1:-1)      = fdr_s(:mkx)  
             cinh_out(i)              = cin_s
             cinlclh_out(i)           = cinlcl_s
             qtten_out(i,mkx:1:-1)    = qtten_s(:mkx)
             slten_out(i,mkx:1:-1)    = slten_s(:mkx)
             ufrc_out(i,mkx:0:-1)     = ufrc_s(0:mkx)
             uflx_out(i,mkx:0:-1)     = uflx_s(0:mkx)  
             vflx_out(i,mkx:0:-1)     = vflx_s(0:mkx)  

             ufrcinvbase_out(i)       = ufrcinvbase_s
             ufrclcl_out(i)           = ufrclcl_s 
             winvbase_out(i)          = winvbase_s
             wlcl_out(i)              = wlcl_s
             plcl_out(i)              = plcl_s
             pinv_out(i)              = pinv_s    
             plfc_out(i)              = plfc_s    
             pbup_out(i)              = pbup_s
             ppen_out(i)              = ppen_s    
             qtsrc_out(i)             = qtsrc_s
             thlsrc_out(i)            = thlsrc_s
             thvlsrc_out(i)           = thvlsrc_s
             emfkbup_out(i)           = emfkbup_s
             cbmflimit_out(i)         = cbmflimit_s
             tkeavg_out(i)            = tkeavg_s
             zinv_out(i)              = zinv_s
             rcwp_out(i)              = rcwp_s
             rlwp_out(i)              = rlwp_s
             riwp_out(i)              = riwp_s

             wu_out(i,mkx:0:-1)       = wu_s(0:mkx)
             qtu_out(i,mkx:0:-1)      = qtu_s(0:mkx)
             thlu_out(i,mkx:0:-1)     = thlu_s(0:mkx)
             thvu_out(i,mkx:0:-1)     = thvu_s(0:mkx)
             uu_out(i,mkx:0:-1)       = uu_s(0:mkx)
             vu_out(i,mkx:0:-1)       = vu_s(0:mkx)
             qtu_emf_out(i,mkx:0:-1)  = qtu_emf_s(0:mkx)
             thlu_emf_out(i,mkx:0:-1) = thlu_emf_s(0:mkx)
             uu_emf_out(i,mkx:0:-1)   = uu_emf_s(0:mkx)
             vu_emf_out(i,mkx:0:-1)   = vu_emf_s(0:mkx)
             uemf_out(i,mkx:0:-1)     = uemf_s(0:mkx)

             dwten_out(i,mkx:1:-1)    = dwten_s(:mkx)
             diten_out(i,mkx:1:-1)    = diten_s(:mkx)
             flxrain_out(i,mkx:0:-1)  = flxrain_s(0:mkx)
             flxsnow_out(i,mkx:0:-1)  = flxsnow_s(0:mkx)
             ntraprd_out(i,mkx:1:-1)  = ntraprd_s(:mkx)
             ntsnprd_out(i,mkx:1:-1)  = ntsnprd_s(:mkx)

             excessu_arr_out(i,mkx:1:-1)  = excessu_arr_s(:mkx)
             excess0_arr_out(i,mkx:1:-1)  = excess0_arr_s(:mkx)
             xc_arr_out(i,mkx:1:-1)       = xc_arr_s(:mkx)
             aquad_arr_out(i,mkx:1:-1)    = aquad_arr_s(:mkx)
             bquad_arr_out(i,mkx:1:-1)    = bquad_arr_s(:mkx)
             cquad_arr_out(i,mkx:1:-1)    = cquad_arr_s(:mkx)
             bogbot_arr_out(i,mkx:1:-1)   = bogbot_arr_s(:mkx)
             bogtop_arr_out(i,mkx:1:-1)   = bogtop_arr_s(:mkx)

             id_exit = .false.
             go to 333

          endif

       endif    

       ! ------------------------------------------------------------------ !
       ! Define a release level, 'prel' and release layer, 'krel'.          !
       ! 'prel' is the lowest level from which buoyancy sorting occurs, and !
       ! 'krel' is the layer index containing 'prel' in it, similar to  the !
       ! previous definitions of 'kinv', 'klcl', and 'klfc'.    In order to !
       ! ensure that only PBL scheme works within the PBL,  if LCL is below !
       ! PBL top height, then 'krel = kinv', while if LCL is above  PBL top !
       ! height, then 'krel = klcl'.   Note however that regardless of  the !
       ! definition of 'krel', cumulus convection induces fluxes within PBL !
       ! through 'fluxbelowinv'.  We can make cumulus convection start from !
       ! any level, even within the PBL by appropriately defining 'krel'  & !
       ! 'prel' here. Then it must be accompanied by appropriate definition !
       ! of source air properties, CIN, and re-setting of 'fluxbelowinv', & !
       ! many other stuffs.                                                 !
       ! Probably the simplest way to adjust 'prel' and 'krel' is to adjust !
       ! 'kinv' because calculation of CIN is based on 'kinv'.              !
       ! Note that even when 'prel' is located above the PBL top height, we !
       ! still have cumulus convection between PBL top height and 'prel':   !
       ! we simply assume that no lateral mixing occurs in this range.      !
       ! ------------------------------------------------------------------ !

       if( klcl.lt.kinv ) then
          krel = kinv
          prel = ps0(krel-1)
          thv0rel = thv0bot(krel) 
       else
          krel = klcl
          prel = plcl 
          thv0rel = thv0lcl
       endif  

       ! --------------------------------------------------------------------------- !
       ! Calculate cumulus base mass flux ('cbmf'), fractional area ('ufrcinv'), and !
       ! and mean vertical velocity (winv) of cumulus updraft at PBL top interface.  !
       ! Also, calculate updraft fractional area (ufrclcl) and vertical velocity  at !
       ! the LCL (wlcl). When LCL is below PBLH, cinlcl = 0 and 'ufrclcl = ufrcinv', !
       ! and 'wlcl = winv.                                                           !
       ! Only updrafts strong enough to overcome CIN can rise over PBL top interface.! 
       ! Thus,  in order to calculate cumulus mass flux at PBL top interface, 'cbmf',!
       ! we need to know 'CIN' ( the strength of potential energy barrier ) and      !
       ! 'sigmaw' ( a standard deviation of updraft vertical velocity at the PBL top !
       ! interface, a measure of turbulentce strength in the PBL ).   Naturally, the !
       ! ratio of these two variables, 'mu' - normalized CIN by TKE- is key variable !
       ! controlling 'cbmf'.  If 'mu' becomes large, only small fraction of updrafts !
       ! with very strong TKE can rise over the PBL - both 'cbmf' and 'ufrc' becomes !
       ! small, but 'winv' becomes large ( this can be easily understood by PDF of w !
       ! at PBL top ).  If 'mu' becomes small, lots of updraft can rise over the PBL !
       ! top - both 'cbmf' and 'ufrc' becomes large, but 'winv' becomes small. Thus, !
       ! all of the key variables associated with cumulus convection  at the PBL top !
       ! - 'cbmf', 'ufrc', 'winv' where 'cbmf = rho*ufrc*winv' - are a unique functi !
       ! ons of 'mu', normalized CIN. Although these are uniquely determined by 'mu',! 
       ! we usually impose two comstraints on 'cbmf' and 'ufrc': (1) because we will !
       ! simply assume that subsidence warming and drying of 'kinv-1' layer in assoc !
       ! iation with 'cbmf' at PBL top interface is confined only in 'kinv-1' layer, !
       ! cbmf must not be larger than the mass within the 'kinv-1' layer. Otherwise, !
       ! instability will occur due to the breaking of CBS condition. If we consider !
       ! semi-Lagrangian vertical advection scheme and explicitly consider the exten !
       ! t of vertical movement of each layer in association with cumulus mass flux, !
       ! we don't need to impose this constraint. However,  using a  semi-Lagrangian !
       ! scheme is a future research subject. Note that this constraint should be ap !
       ! plied for all interfaces above PBL top as well as PBL top interface.   As a !
       ! result, this 'cbmf' constraint impose a 'lower' limit on mu - 'mumin0'. (2) !
       ! in order for mass flux parameterization - rho*(w'a')= M*(a_c-a_e) - to   be !
       ! valid, cumulus updraft fractional area should be much smaller than 1.    In !
       ! current code, we impose 'rmaxfrac = 0.1 ~ 0.2'   through the whole vertical !
       ! layers where cumulus convection occurs. At the PBL top interface,  the same !
       ! constraint is made by imposing another lower 'lower' limit on mu, 'mumin1'. !
       ! After that, also limit 'ufrclcl' to be smaller than 'rmaxfrac' by 'mumin2'. !
       ! --------------------------------------------------------------------------- !
       
       ! --------------------------------------------------------------------------- !
       ! Calculate normalized CIN, 'mu' satisfying all the three constraints imposed !
       ! on 'cbmf'('mumin0'), 'ufrc' at the PBL top - 'ufrcinv' - ( by 'mumin1' from !
       ! a parameter sentence), and 'ufrc' at the LCL - 'ufrclcl' ( by 'mumin2').    !
       ! Note that 'cbmf' does not change between PBL top and LCL  because we assume !
       ! that buoyancy sorting does not occur when cumulus updraft is unsaturated.   !
       ! --------------------------------------------------------------------------- !
   
       if(use_CINcin) then       
          wcrit   = sqrt( 2._r8 * cin * rbuoy )      
       else
          wcrit   = sqrt( 2._r8 * cinlcl * rbuoy )   
       endif
       sigmaw  = sqrt( rkfre * tkeavg + epsvarw )
       mu  = wcrit/sigmaw/1.4142_r8                  
       if(mu .ge. 3._r8) then
        ! write(iulog,*) 'mu >= 3'
          id_exit = .true.
          go to 333
       endif
       rho0inv = ps0(kinv-1)/(r*thv0top(kinv-1)*exns0(kinv-1))
       cbmf = (rho0inv*sigmaw/2.5066_r8)*exp(-mu**2)
       ! 1. 'cbmf' constraint
       cbmflimit = 0.9_r8*dp0(kinv-1)/g/dt
       mumin0 = 0._r8
       if(cbmf.gt.cbmflimit) mumin0 = sqrt(-log(2.5066_r8*cbmflimit/rho0inv/sigmaw))
       ! 2. 'ufrcinv' constraint
       mu = max(max(mu,mumin0),mumin1)
       ! 3. 'ufrclcl' constraint      
       mulcl = sqrt(2.*cinlcl*rbuoy)/1.4142/sigmaw
       mulclstar = sqrt(max(0._r8,2*(exp(-mu**2)/2.5066_r8)**2*(1/erfc(mu)**2-0.25_r8/rmaxfrac**2)))
       if(mulcl.gt.1.e-8_r8.and.mulcl.gt.mulclstar) then
          mumin2 = compute_mumin2(mulcl,rmaxfrac,mu)
          if(mu.gt.mumin2) then
             write(iulog,*) 'Critical error in mu calculation in UW_ShCu'
             stop
          endif
          mu = max(mu,mumin2)
          if(mu.eq.mumin2) limit_ufrc(i) = 1._r8
       endif
       if(mu.eq.mumin0) limit_cbmf(i) = 1._r8
       if(mu.eq.mumin1) limit_ufrc(i) = 1._r8

       ! ------------------------------------------------------------------- !    
       ! Calculate final ['cbmf','ufrcinv','winv'] at the PBL top interface. !
       ! Note that final 'cbmf' here is obtained in such that 'ufrcinv' and  !
       ! 'ufrclcl' are smaller than ufrcmax and there is no CBF instability. !
       ! ------------------------------------------------------------------- !

       cbmf = (rho0inv*sigmaw/2.5066_r8)*exp(-mu**2)                       
       winv = sigmaw*(2/2.5066_r8)*exp(-mu**2)/erfc(mu)
       ufrcinv = cbmf/winv/rho0inv

       ! ------------------------------------------------------------------- !
       ! Calculate ['ufrclcl','wlcl'] at the LCL. When LCL is below PBL top, !
       ! it automatically becomes 'ufrclcl = ufrcinv' & 'wlcl = winv', since !
       ! it was already set to 'cinlcl=0' if LCL is below PBL top interface. !
       ! Note 'cbmf' at the PBL top is the same as 'cbmf' at the LCL.  Note  !
       ! also that final 'cbmf' here is obtained in such that 'ufrcinv' and  !
       ! 'ufrclcl' are smaller than ufrcmax and there is no  CBF instability.!
       ! By construction, it must be 'wlcl2>0' but for assurance, I checked  !
       ! this again in the below block. If 'ufrclcl < 0.1%', just exit.      !
       ! ------------------------------------------------------------------- !

       wtw = winv * winv - 2._r8 * cinlcl * rbuoy
       if(wtw .le. 0._r8) then
        ! write(iulog,*) 'wlcl2 < 0 at the LCL'
          exit_wtw(i) = 1._r8
          id_exit = .true.
          go to 333
       endif
       wlcl = sqrt(wtw)
       ufrclcl = cbmf/wlcl/rho0inv
       wrel = wlcl
       if(ufrclcl.le.0.001_r8) then
        ! write(iulog,*) 'ufrclcl <= 0.001' 
          exit_ufrc(i) = 1._r8
          id_exit = .true.
          go to 333
       endif
       ufrc(krel-1) = ufrclcl

       ! ----------------------------------------------------------------------- !
       ! Below is just diagnostic output for detailed analysis of cumulus scheme !
       ! ----------------------------------------------------------------------- !

       ufrcinvbase = ufrcinv
       winvbase = winv
       umf(kinv-1:krel-1) = cbmf   
       wu(kinv-1:krel-1) = winv   

       ! -------------------------------------------------------------------------- ! 
       ! Define updraft properties at the level where buoyancy sorting starts to be !
       ! happening, i.e., by definition, at 'prel' level within the release layer.  !
       ! Because no lateral entrainment occurs upto 'prel', conservative scalars of ! 
       ! cumulus updraft at release level is same as those of source air.  However, ! 
       ! horizontal momentums of source air are modified by horizontal PGF forcings ! 
       ! from PBL top interface to 'prel'.  For this case, we should add additional !
       ! horizontal momentum from PBL top interface to 'prel' as will be done below !
       ! to 'usrc' and 'vsrc'. Note that below cumulus updraft properties - umf, wu,!
       ! thlu, qtu, thvu, uu, vu - are defined all interfaces not at the layer mid- !
       ! point. From the index notation of cumulus scheme, wu(k) is the cumulus up- !
       ! draft vertical velocity at the top interface of k layer.                   !
       ! -------------------------------------------------------------------------- !

       emf(krel-1)  = 0._r8
       umf(krel-1)  = cbmf
       wu(krel-1)   = wrel
       thlu(krel-1) = thlsrc
       qtu(krel-1)  = qtsrc
       call conden(prel,thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check,qsat)
       if(id_check.eq.1) then
          exit_conden(i) = 1._r8
          id_exit = .true.
          go to 333
       end if
       thvu(krel-1) = thj*(1._r8 + zvir*qvj - qlj - qij )       

       uplus = 0._r8
       vplus = 0._r8
       if( krel.eq.kinv ) then
           uplus = PGFc * ssu0(kinv) * ( prel - ps0(kinv-1) )
           vplus = PGFc * ssv0(kinv) * ( prel - ps0(kinv-1) )
       else
         do k = kinv, max(krel-1,kinv)
           uplus = uplus + PGFc * ssu0(k) * ( ps0(k) - ps0(k-1) )
           vplus = vplus + PGFc * ssv0(k) * ( ps0(k) - ps0(k-1) )
         end do
           uplus = uplus + PGFc * ssu0(krel) * ( prel - ps0(krel-1) )
           vplus = vplus + PGFc * ssv0(krel) * ( prel - ps0(krel-1) )
       end if
       uu(krel-1) = usrc + uplus
       vu(krel-1) = vsrc + vplus      

       ! -------------------------------------------------------------------------- !
       ! Define environmental properties at the level where buoyancy sorting occurs !
       ! ('pe', normally, layer midpoint except in the 'krel' layer). In the 'krel' !
       ! layer where buoyancy sorting starts to occur, however, 'pe' is defined     !
       ! differently because LCL is regarded as lower interface for mixing purpose. !
       ! -------------------------------------------------------------------------- !

       pe      = 0.5_r8 * ( prel + ps0(krel) )
       dpe     = prel - ps0(krel)
       exne    = exnf(pe)
       thvebot = thv0rel
       thle    = thl0(krel) + ssthl0(krel) * ( pe - p0(krel) )
       qte     = qt0(krel)  + ssqt0(krel)  * ( pe - p0(krel) )
       ue      = u0(krel)   + ssu0(krel)   * ( pe - p0(krel) )
       ve      = v0(krel)   + ssv0(krel)   * ( pe - p0(krel) )

       !-------------------------! 
       ! Buoyancy-Sorting Mixing !
       !-------------------------!------------------------------------------------ !
       !                                                                           !
       !  In order to complete buoyancy-sorting mixing at layer mid-point, and so  ! 
       !  calculate 'updraft mass flux, updraft w velocity, conservative scalars'  !
       !  at the upper interface of each layer, we need following 3 information.   ! 
       !                                                                           !
       !  1. Pressure where mixing occurs ('pe'), and temperature at 'pe' which is !
       !     necessary to calculate various thermodynamic coefficients at pe. This !
       !     temperature is obtained by undiluted cumulus properties lifted to pe. ! 
       !  2. Undiluted updraft properties at pe - conservative scalar and vertical !
       !     velocity -which are assumed to be the same as the properties at lower !
       !     interface only for calculation of fractional lateral entrainment  and !
       !     detrainment rate ( fer(k) and fdr(k) [Pa-1] ), respectively.    Final !
       !     values of cumulus conservative scalars and w at the top interface are !
       !     calculated afterward after obtaining fer(k) & fdr(k).                 !
       !  3. Environmental properties at pe.                                       !
       ! ------------------------------------------------------------------------- !
       
       ! ------------------------------------------------------------------------ ! 
       ! Define cumulus scale height.                                             !
       ! Cumulus scale height is defined as the maximum height cumulus can reach. !
       ! In case of premitive code, cumulus scale height ('cush')  at the current !
       ! time step was assumed to be the same as 'cush' of previous time step.    !
       ! However, I directly calculated cush at each time step using an iterative !
       ! method. Note that within the cumulus scheme, 'cush' information is  used !
       ! only at two places during buoyancy-sorting process:                      !
       ! (1) Even negatively buoyancy mixtures with strong vertical velocity      !
       !     enough to rise up to 'rle*scaleh' (rle = 0.1) from pe are entrained  !
       !     into cumulus updraft,                                                !  
       ! (2) The amount of mass that is involved in buoyancy-sorting mixing       !
       !      process at pe is rei(k) = rkm/scaleh/rho*g [Pa-1]                   !
       ! In terms of (1), I think critical stopping distance might be replaced by !
       ! layer thickness. In future, we will use rei(k) = (0.5*rkm/z0(k)/rho/g).  !
       ! In the premitive code,  'scaleh' was largely responsible for the jumping !
       ! variation of precipitation amount.                                       !
       ! ------------------------------------------------------------------------ !   

       scaleh = tscaleh
       if(tscaleh .lt. 0.0_r8) scaleh = 1000._r8 

       do iter_scaleh = 1, 3

       ! ---------------------------------------------------------------- !
       ! Initialization of 'kbup' and 'kpen'                              !
       ! ---------------------------------------------------------------- !
       ! 'kbup' is the top-most layer in which cloud buoyancy is positive !
       ! both at the top and bottom interface of the layer. 'kpen' is the !
       ! layer upto which cumulus panetrates ,i.e., cumulus w at the base !
       ! interface is positive, but becomes negative at the top interface.!
       ! Here, we initialize 'kbup' and 'kpen'. These initializations are !  
       ! not trivial but important, expecially   in calculating turbulent !
       ! fluxes without confliction among several physics as explained in !
       ! detail in the part of turbulent fluxes calculation later.   Note !
       ! that regardless of whether 'kbup' and 'kpen' are updated or  not !
       ! during updraft motion,  penetrative entrainments are dumped down !
       ! across the top interface of 'kbup' later.      More specifically,!
       ! penetrative entrainment heat and moisture fluxes are  calculated !
       ! from the top interface of 'kbup' layer  to the base interface of !
       ! 'kpen' layer. Because of this, initialization of 'kbup' & 'kpen' !
       ! influence the convection system when there are not updated.  The !  
       ! below initialization of 'kbup = krel' assures  that  penetrative !
       ! entrainment fluxes always occur at interfaces above the PBL  top !
       ! interfaces (i.e., only at interfaces k >=kinv ), which seems  to !
       ! be attractable considering that the most correct fluxes  at  the !
       ! PBL top interface can be ontained from the 'fluxbelowinv'  using !
       ! reconstructed PBL height.                                        ! 
       ! The 'kbup = krel'(after going through the whole buoyancy sorting !
       ! proces during updraft motion) implies that cumulus updraft  from !
       ! the PBL top interface can not reach to the LFC,so that 'kbup' is !
       ! not updated during upward. This means that cumulus updraft   did !
       ! not fully overcome the buoyancy barrier above just the PBL top.  !
       ! If 'kpen' is not updated either ( i.e., cumulus cannot rise over !
       ! the top interface of release layer),penetrative entrainment will !
       ! not happen at any interfaces.  If cumulus updraft can rise above !
       ! the release layer but cannot fully overcome the buoyancy barrier !
       ! just above PBL top interface, penetratve entrainment   occurs at !
       ! several above interfaces, including the top interface of release ! 
       ! layer. In the latter case, warming and drying tendencies will be !
       ! be initiated in 'krel' layer. Note current choice of 'kbup=krel' !
       ! is completely compatible with other flux physics without  double !
       ! or miss counting turbulent fluxes at any interface. However, the !
       ! alternative choice of 'kbup=krel-1' also has itw own advantage - !
       ! when cumulus updraft cannot overcome buoyancy barrier just above !
       ! PBL top, entrainment warming and drying are concentrated in  the !
       ! 'kinv-1' layer instead of 'kinv' layer for this case. This might !
       ! seems to be more dynamically reasonable, but I will choose the   !
       ! 'kbup = krel' choice since it is more compatible  with the other !
       ! parts of the code, expecially, when we chose ' use_emf=.false. ' !
       ! as explained in detail in turbulent flux calculation part.       !
       ! ---------------------------------------------------------------- ! 

       kbup    = krel
       kpen    = krel
       
       ! ------------------------------------------------------------ !
       ! Since 'wtw' is continuously updated during vertical motion,  !
       ! I need below initialization command within this 'iter_scaleh'!
       ! do loop. Similarily, I need initializations of environmental !
       ! properties at 'krel' layer as below.                         !
       ! ------------------------------------------------------------ !

       wtw     = wlcl * wlcl
       pe      = 0.5_r8 * ( prel + ps0(krel) )
       dpe     = prel - ps0(krel)
       exne    = exnf(pe)
       thvebot = thv0rel
       thle    = thl0(krel) + ssthl0(krel) * ( pe - p0(krel) )
       qte     = qt0(krel)  + ssqt0(krel)  * ( pe - p0(krel) )
       ue      = u0(krel)   + ssu0(krel)   * ( pe - p0(krel) )
       ve      = v0(krel)   + ssv0(krel)   * ( pe - p0(krel) )

       ! ----------------------------------------------------------------------- !
       ! Cumulus rises upward from 'prel' ( or base interface of  'krel' layer ) !
       ! until updraft vertical velocity becomes zero.                           !
       ! Buoyancy sorting is performed via two stages. (1) Using cumulus updraft !
       ! properties at the base interface of each layer,perform buoyancy sorting !
       ! at the layer mid-point, 'pe',  and update cumulus properties at the top !
       ! interface, and then  (2) by averaging updated cumulus properties at the !
       ! top interface and cumulus properties at the base interface,   calculate !
       ! cumulus updraft properties at pe that will be used  in buoyancy sorting !
       ! mixing - thlue, qtue and, wue.  Using this averaged properties, perform !
       ! buoyancy sorting again at pe, and re-calculate fer(k) and fdr(k). Using !
       ! this recalculated fer(k) and fdr(k),  finally calculate cumulus updraft !
       ! properties at the top interface - thlu, qtu, thvu, uu, vu. In the below,!
       ! 'iter_xc = 1' performs the first stage, while 'iter_xc= 2' performs the !
       ! second stage. We can increase the number of iterations, 'nter_xc'.as we !
       ! want, but a sample test indicated that about 3 - 5 iterations  produced !
       ! satisfactory converent solution. Finally, identify 'kbup' and 'kpen'.   !
       ! ----------------------------------------------------------------------- !
       
       do k = krel, mkx - 1 ! Here, 'k' is a layer index.

          km1 = k - 1

          thlue = thlu(km1)
          qtue  = qtu(km1)    
          wue   = wu(km1)
          wtwb  = wtw  

        do iter_xc = 1, niter_xc
          
          wtw = wu(km1) * wu(km1)

          ! ---------------------------------------------------------------- !
          ! Calculate environmental and cumulus saturation 'excess' at 'pe'. !
          ! Note that in order to calculate saturation excess, we should use ! 
          ! liquid water temperature instead of temperature  as the argument !
          ! of "qsat". However, normal argument of "qsat" is temperature.    !
          ! ---------------------------------------------------------------- !

          call conden(pe,thle,qte,thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
          end if
          thv0j    = thj * ( 1._r8 + zvir*qvj - qlj - qij )
          rho0j    = pe / ( r * thv0j * exne )
          qsat_arg = thle*exne     
          status   = qsat(qsat_arg,pe,es(1),qs(1),gam(1),1)
          excess0  = qte - qs(1)

          call conden(pe,thlue,qtue,thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
          end if
          ! ----------------------------------------------------------------- !
          ! Detrain excessive condensate larger than 'criqc' from the cumulus ! 
          ! updraft before performing buoyancy sorting. All I should to do is !
          ! to update 'thlue' &  'que' here. Below modification is completely !
          ! compatible with the other part of the code since 'thule' & 'qtue' !
          ! are used only for buoyancy sorting. I found that as long as I use !
          ! 'niter_xc >= 2',  detraining excessive condensate before buoyancy !
          ! sorting has negligible influence on the buoyancy sorting results. !   
          ! ----------------------------------------------------------------- !
          if( (qlj + qij) .gt. criqc ) then
             exql  = ( ( qlj + qij ) - criqc ) * qlj / ( qlj + qij )
             exqi  = ( ( qlj + qij ) - criqc ) * qij / ( qlj + qij )
             qtue  = qtue - exql - exqi
             thlue = thlue + (xlv/cp/exne)*exql + (xls/cp/exne)*exqi 
          endif
          call conden(pe,thlue,qtue,thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
          end if
          thvj     = thj * ( 1._r8 + zvir * qvj - qlj - qij )
          tj       = thj * exne ! This 'tj' is used for computing thermo. coeffs. below
          qsat_arg = thlue*exne
          status   = qsat(qsat_arg,pe,es(1),qs(1),gam(1),1)
          excessu  = qtue - qs(1)

          ! ------------------------------------------------------------------- !
          ! Calculate critical mixing fraction, 'xc'. Mixture with mixing ratio !
          ! smaller than 'xc' will be entrained into cumulus updraft.  Both the !
          ! saturated updrafts with 'positive buoyancy' or 'negative buoyancy + ! 
          ! strong vertical velocity enough to rise certain threshold distance' !
          ! are kept into the updraft in the below program. If the core updraft !
          ! is unsaturated, we can set 'xc = 0' and let the cumulus  convection !
          ! still works or we may exit.                                         !
          ! Current below code does not entrain unsaturated mixture. However it !
          ! should be modified such that it also entrain unsaturated mixture.   !
          ! ------------------------------------------------------------------- !

          ! ----------------------------------------------------------------- !
          ! cridis : Critical stopping distance for buoyancy sorting purpose. !
          ! ----------------------------------------------------------------- !

            cridis = rle*scaleh              ! Original code
          ! cridis = 1._r8*(zs0(k) - zs0(k-1))  ! New code
 
       if(use_newbuosort) then

          ! ------------------------------ !
          ! New Buoyancy Sorting Algorithm !
          ! ------------------------------ !                   

          ! ----------------------------------------------------------------- !
          ! Case 1 : When both cumulus and env. are unsaturated or saturated. !
          ! ----------------------------------------------------------------- !
          if((excessu.lt.0.and.excess0.lt.0).or.(excessu.gt.0.and.excess0.gt.0)) then
              xc = min(1._r8,max(0._r8,1._r8-2._r8*rbuoy*g*cridis/wue**2._r8*(1-thvj/thv0j)))
              ! Below 3 lines are diagnostic output not influencing
              ! numerical calculations.
              aquad = 0._r8
              bquad = 0._r8
              cquad = 0._r8
          else
          ! -------------------------------------------------- !
          ! Case 2 : When either cumulus or env. is saturated. !
          ! -------------------------------------------------- !
              xsat = excessu / ( excessu - excess0 );
              thlxsat = thlue + xsat * ( thle - thlue );
              qtxsat = qtue + xsat * ( qte - qtue );
              call conden(pe,thlxsat,qtxsat,thj,qvj,qlj,qij,qse,id_check,qsat)
              if(id_check.eq.1) then
                 exit_conden(i) = 1._r8
                 id_exit = .true.
                 go to 333
              end if
              thvxsat = thj * ( 1._r8 + zvir * qvj - qlj - qij )               
              ! -------------------------------------------------- !
              ! kk=1 : Cumulus Segment, kk=2 : Environment Segment !
              ! -------------------------------------------------- ! 
              do kk = 1, 2 
                   if( kk.eq.1 ) then
                       thv_x0 = thvj;
                       thv_x1 = ( 1._r8 - 1._r8/xsat ) * thvj + ( 1._r8/xsat ) * thvxsat;
                   else
                       thv_x1 = thv0j;
                       thv_x0 = ( xsat / ( xsat - 1._r8 ) ) * thv0j + ( 1._r8/( 1._r8 - xsat ) ) * thvxsat;
                   endif
                   aquad =  wue**2;
                   bquad =  2._r8*rbuoy*g*cridis*(thv_x1 - thv_x0)/thv0j - 2._r8*wue**2;
                   cquad =  2._r8*rbuoy*g*cridis*(thv_x0 - thv0j) /thv0j +    wue**2;
                   if( kk.eq.1 ) then
                       if((bquad**2-4._r8*aquad*cquad).ge.0._r8) then
                           call roots(aquad,bquad,cquad,xs1,xs2,status)
                           x_cu = min(1._r8,max(0._r8,min(xsat,min(xs1,xs2))))
                       else
                           x_cu = xsat;
                       endif
                   else 
                       if((bquad**2-4._r8*aquad*cquad).ge.0._r8) then
                           call roots(aquad,bquad,cquad,xs1,xs2,status)
                           x_en = min(1._r8,max(0._r8,max(xsat,min(xs1,xs2))))
                       else
                           x_en = 1._r8;
                       endif
                   endif
              enddo
              if( x_cu.eq.xsat ) then
                 xc = max(x_cu, x_en);
              else
                 xc = x_cu;
              endif
          endif

       else
 
          ! ------------------------------ !
          ! Old Buoyancy Sorting Algorithm !
          ! ------------------------------ !

          if( excessu.lt.0 ) then                       
              ! ----------------------------------------------------- !
              ! Ideally, it should be 'xc = 1'   when 'thvj > thv0j', !
              !                       'xc = 0'   when 'thvj < thv0j', !
              !                       'xc = 0.5' when 'thvj = thv0j'  !
              ! This should be modified layer. Also, 'exit_drycore'   !
              ! should be removed.                                    !
              ! ----------------------------------------------------- !
              xc = 0.
            ! write(iulog,*) 'Dry Core Updraft Warning in UW-shallow scheme'
              exit_drycore(i) = 1._r8
            ! id_exit = .true.
            ! go to 333 
          elseif( excessu.ge.0.and.excess0.ge.0 ) then
              ! ---------------------------------------------------------------- !
              ! When both cumulus and environment are saturated,  I should check !
              ! if below calculation assuming that resulting buoyancy of mixture !
              ! is a simple linear function of 'x' is correctly or not.          !
              ! ---------------------------------------------------------------- !
              aquad =  wue**2
              bquad = -2._r8*rbuoy*g*cridis*(thvj - thv0j)/thv0j - 2._r8*wue**2
              cquad =  2._r8*rbuoy*g*cridis*(thvj - thv0j)/thv0j +    wue**2
              if( (-2._r8*aquad-bquad).ge.0._r8 ) then
                xc = 1._r8
              else
                xc = min(1._r8,max(0._r8,-bquad/aquad-1._r8)) 
              endif
          elseif ( excessu.ge.0.and.excess0.lt.0 ) then
              ! --------------------------------------------------------------------- !
              ! Below code does not entrain 'unsaturated positively buoyant mixture', !
              ! which is clearly wrong. I should modify current below code such that  !
              ! cumulus also entrain 'unsaturated positively buoyant mixture'.        !
              ! --------------------------------------------------------------------- !
              xsat   = excessu / ( excessu - excess0 )
              rdqsdt = ep2 * xlv * qse / r / tj**2
              rbeta  = ( 1._r8 + ( 1._r8 + zvir ) * tj * rdqsdt ) / ( 1._r8 + rdqsdt * xlv / cp )
              ths0   = min(thv0j,thvj + rbeta*(thle - thlue) + (rbeta*xlv/cp/exne - thj)*(qte - qtue))
              aquad =  wue**2
              bquad = -2._r8*rbuoy*g*cridis*(thvj - ths0) /thv0j - 2._r8*wue**2
              cquad =  2._r8*rbuoy*g*cridis*(thvj - thv0j)/thv0j +    wue**2
              if((bquad**2-4._r8*aquad*cquad).ge.0._r8) then
                 ! ------------------------------------------------------------------- !
                 ! I checked that below premitive code produced almost ( not exactly ) !
                 ! same results as mine. For numerical stability in solving quadrature !
                 ! equation let's use the primitive.                                   !
                 ! ------------------------------------------------------------------- ! 
                 call roots(aquad,bquad,cquad,xs1,xs2,status) ! Premitive
                 xc = min(1._r8,max(0._r8,min(xsat,min(xs1,xs2))))  ! Premitive
               ! xc = min(1._r8,max(0._r8,min(xsat,(-bquad-sqrt(bquad**2-4._r8*aquad*cquad))/(2._r8*aquad)))) ! Sungsu
              else               
                 xc = xsat
              endif
          endif

       endif

          ! ------------------------------------------------------------------------ !
          ! Compute fractional lateral entrainment & detrainment rate in each layers.!
          ! The unit of rei(k), fer(k), and fdr(k) is [Pa-1].  Alternative choice of !
          ! 'rei(k)' is also shown below, where coefficient 0.5 was from approximate !
          ! tuning against the BOMEX case.                                           !
          ! In order to prevent the onset of instability in association with cumulus !
          ! induced subsidence advection, cumulus mass flux at the top interface  in !
          ! any layer should be smaller than ( 90% of ) total mass within that layer.!
          ! I imposed limits on 'rei(k)' as below,  in such that stability condition ! 
          ! is always satisfied.                                                     !
          ! Below limiter of 'rei(k)' becomes negative for some cases, causing error.!
          ! So, for the time being, I came back to the original limiter.             !
          ! ------------------------------------------------------------------------ !

          ee2    = xc**2
          ud2    = 1._r8 - 2.*xc + xc**2
        ! rei(k) = (rkm/scaleh/g/rho0j)    ! Default.
          rei(k) = (0.5_r8*rkm/z0(k)/g/rho0j) ! Alternative.
          ! ------------------------------ !
          ! Below is the original limiter. !
          ! ------------------------------ !
          if( xc.gt.0.5_r8 ) rei(k) = min(rei(k),0.9*log(dp0(k)/g/dt/umf(km1) + 1._r8)/dpe/(2._r8*xc-1._r8))
          ! ------------------------------ !
          ! Below is the new limiter.      !
          ! ------------------------------ !
          ! if(xc.gt.0.5_r8) then
          !    limit_rei(i) = 1._r8
          !    rei(k) = min(rei(k),log(0.9*dp0(k)/g/dt/umf(km1))/dpe/(2._r8*xc-1._r8))
          ! elseif(xc.eq.0.5_r8) then
          !    if((0.9_r8*dp0(k)/g/dt/umf(km1)).lt.1._r8) then
          !        write(iulog,*) 'Impossible to avoid CBS instability in mcshallow.F90'
          !        id_exit = .true.
          !        go to 333           
          !        exit_rei(i) = 1._r8
          !    endif
          ! else
          !    if((0.9*dp0(k)/g/dt/umf(km1)).lt.1._r8) then 
          !        limit_rei(i) = 1._r8
          !        rei(k) = max(rei(k),log(0.9_r8*dp0(k)/g/dt/umf(km1))/dpe/(2._r8*xc-1.))
          !    endif
          ! endif
          fer(k) = rei(k) * ee2
          fdr(k) = rei(k) * ud2

          ! ------------------------------------------------------------------------------ !
          ! Iteration Start due to 'maxufrc' constraint [ ****************************** ] ! 
          ! ------------------------------------------------------------------------------ !

          ! -------------------------------------------------------------------------- !
          ! Calculate cumulus updraft mass flux and penetrative entrainment mass flux. !
          ! Note that  non-zero penetrative entrainment mass flux will be asigned only !
          ! to interfaces from the top interface of 'kbup' layer to the base interface !
          ! of 'kpen' layer as will be shown later.                                    !
          ! -------------------------------------------------------------------------- !

          umf(k) = umf(km1) * exp( dpe * ( fer(k) - fdr(k) ) )
          emf(k) = 0.0_r8    

          ! --------------------------------------------------------- !
          ! Compute cumulus updraft properties at the top interface.  !
          ! Also use Tayler expansion in order to treat limiting case !
          ! --------------------------------------------------------- !

          if( fer(k)*dpe .lt. 1.e-4_r8 ) then
            thlu(k) = thlu(km1) + ( thle + ssthl0(k) * dpe / 2._r8 - thlu(km1) ) * fer(k) * dpe
            qtu(k)  =  qtu(km1) + ( qte  + ssqt0(k) * dpe / 2._r8  -  qtu(km1) ) * fer(k) * dpe
            uu(k)   =   uu(km1) + ( ue   + ssu0(k) * dpe / 2._r8   -   uu(km1) ) * fer(k) * dpe - PGFc * ssu0(k) * dpe
            vu(k)   =   vu(km1) + ( ve   + ssv0(k) * dpe / 2._r8   -   vu(km1) ) * fer(k) * dpe - PGFc * ssv0(k) * dpe
          else
            thlu(k) = ( thle + ssthl0(k) / fer(k) - ssthl0(k) * dpe / 2._r8 ) -          &
                      ( thle + ssthl0(k) * dpe / 2._r8 - thlu(km1) + ssthl0(k) / fer(k) ) * exp(-fer(k) * dpe)
            qtu(k)  = ( qte  + ssqt0(k) / fer(k) - ssqt0(k) * dpe / 2._r8 ) -            &  
                      ( qte  + ssqt0(k) * dpe / 2._r8 - qtu(km1) + ssqt0(k) / fer(k) ) * exp(-fer(k) * dpe)
            uu(k) =   ( ue + ( 1 - PGFc ) * ssu0(k) / fer(k) - ssu0(k) * dpe / 2._r8 ) - &
                      ( ue + ssu0(k) * dpe / 2._r8 - uu(km1) + ( 1 - PGFc ) * ssu0(k) / fer(k) ) * exp(-fer(k) * dpe)
            vu(k) =   ( ve + ( 1 - PGFc ) * ssv0(k) / fer(k) - ssv0(k) * dpe / 2._r8 ) - &
                      ( ve + ssv0(k) * dpe / 2._r8 - vu(km1) + ( 1 - PGFc ) * ssv0(k) / fer(k) ) * exp(-fer(k) * dpe)
          end if

          !------------------------------------------------------------------- !
          ! Expel some of cloud water and ice from cumulus  updraft at the top !
          ! interface.  Note that this is not 'detrainment' term  but a 'sink' !
          ! term of cumulus updraft qt ( or one part of 'source' term of  mean !
          ! environmental qt ). At this stage, as the most simplest choice, if !
          ! condensate amount within cumulus updraft is larger than a critical !
          ! value, 'criqc', expels the surplus condensate from cumulus updraft !
          ! to the environment. A certain fraction ( e.g., 'frc_sus' ) of this !
          ! expelled condesnate will be in a form that can be suspended in the !
          ! layer k where it was formed, while the other fraction, '1-frc_sus' ! 
          ! will be in a form of precipitatble (e.g.,can potentially fall down !
          ! across the base interface of layer k ). In turn we should describe !
          ! subsequent falling of precipitable condensate ('1-frc_sus') across !
          ! the base interface of the layer k, &  evaporation of precipitating !
          ! water in the below layer k-1 and associated evaporative cooling of !
          ! the later, k-1, and falling of 'non-evaporated precipitating water !
          ! ( which was initially formed in layer k ) and a newly-formed preci !
          ! pitable water in the layer, k-1', across the base interface of the !
          ! lower layer k-1.  Cloud microphysics should correctly describe all !
          ! of these process.  In a near future, I should significantly modify !
          ! this cloud microhpysics, including precipitation-induced downdraft !
          ! also.                                                              !
          ! ------------------------------------------------------------------ !

          call conden(ps0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
          end if
          if( (qlj + qij) .gt. criqc ) then
             exql    = ( ( qlj + qij ) - criqc ) * qlj / ( qlj + qij )
             exqi    = ( ( qlj + qij ) - criqc ) * qij / ( qlj + qij )
             ! ---------------------------------------------------------------- !
             ! It is very important to re-update 'qtu' and 'thlu'  at the upper ! 
             ! interface after expelling condensate from cumulus updraft at the !
             ! top interface of the layer. As mentioned above, this is a 'sink' !
             ! of cumulus qt (or equivalently, a 'source' of environmentasl qt),!
             ! not a regular convective'detrainment'.                           !
             ! ---------------------------------------------------------------- !
             qtu(k)  = qtu(k) - exql - exqi
             thlu(k) = thlu(k) + (xlv/cp/exns0(k))*exql + (xls/cp/exns0(k))*exqi 
             ! ---------------------------------------------------------------- !
             ! Expelled cloud condensate into the environment from the updraft. ! 
             ! After all the calculation later, 'dwten' and 'diten' will have a !
             ! unit of [ kg/kg/s ], because it is a tendency of qt. Restoration !
             ! of 'dwten' and 'diten' to this correct unit through  multiplying !
             ! 'umf(k)*g/dp0(k)' will be performed later after finally updating !
             ! 'umf' using a 'rmaxfrac' constraint near the end of this updraft !
             ! buoyancy sorting loop.                                           !
             ! ---------------------------------------------------------------- !
             dwten(k) = exql   
             diten(k) = exqi
          else
             dwten(k) = 0._r8
             diten(k) = 0._r8
          endif
          ! ----------------------------------------------------------------- ! 
          ! Update 'thvu(k)' after detraining condensate from cumulus updraft.!
          ! ----------------------------------------------------------------- ! 
          call conden(ps0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
          end if  
          thvu(k) = thj * ( 1._r8 + zvir * qvj - qlj - qij )

          ! ----------------------------------------------------------- ! 
          ! Calculate updraft vertical velocity at the upper interface. !
          ! In order to calculate 'wtw' at the upper interface, we use  !
          ! 'wtw' at the lower interface. Note  'wtw'  is continuously  ! 
          ! updated as cumulus updraft rises.                           !
          ! ----------------------------------------------------------- !

          bogbot = rbuoy * ( thvu(km1) / thvebot  - 1._r8 ) ! Cloud buoyancy at base interface
          bogtop = rbuoy * ( thvu(k) / thv0top(k) - 1._r8 ) ! Cloud buoyancy at top  interface

          delbog = bogtop - bogbot
          drage  = fer(k) * ( 1._r8 + rdrag )
          expfac = exp(-2._r8*drage*dpe)

          wtwb = wtw
          if( drage*dpe .gt. 1.e-3_r8 ) then
             wtw = wtw*expfac + (delbog + (1._r8-expfac)*(bogbot + delbog/(-2._r8*drage*dpe)))/(rho0j*drage)
          else
             wtw = wtw + dpe * ( bogbot + bogtop ) / rho0j
          endif
         
          ! -------------------------------------------------------------- !
          ! Repeat 'iter_xc' iteration loop until 'iter_xc = niter_xc'.    !
          ! Also theat the case even when wtw < 0 at the 'kpen' interface. !
          ! -------------------------------------------------------------- !  
          
          if(wtw.gt.0._r8) then   
             thlue = 0.5_r8 * (thlu(km1)+thlu(k))
             qtue  = 0.5_r8 * (qtu(km1)+qtu(k))         
             wue   = 0.5_r8 * sqrt(max(wtwb+wtw,0._r8))
          else
             go to 111
          endif 

       enddo ! End of 'iter_xc' loop  

   111 continue

          ! --------------------------------------------------------------------------- ! 
          ! Add the contribution of self-detrainment  to vertical variations of cumulus !
          ! updraft mass flux. The reason why we are trying to include self-detrainment !
          ! is as follows.  In current scheme,  vertical variation of updraft mass flux !
          ! is not fully consistent with the vertical variation of updraft vertical w.  !
          ! For example, within a given layer, let's assume that  cumulus w is positive !
          ! at the base interface, while negative at the top interface. This means that !
          ! cumulus updraft cannot reach to the top interface of the layer. However,    !
          ! cumulus updraft mass flux at the top interface is not zero according to the !
          ! vertical tendency equation of cumulus mass flux.   Ideally, cumulus updraft ! 
          ! mass flux at the top interface should be zero for this case. In order to    !
          ! assures that cumulus updraft mass flux goes to zero when cumulus updraft    ! 
          ! vertical velocity goes to zero, we are imposing self-detrainment term as    !
          ! below by considering layer-mean cloud buoyancy and cumulus updraft vertical !
          ! velocity square at the top interface. Use of auto-detrainment term will  be !
          ! determined by setting 'use_self_detrain=.true.' in the parameter sentence.  !
          ! --------------------------------------------------------------------------- !
     
          if(use_self_detrain) then
             autodet = min( 0.5_r8*g*(bogbot+bogtop)/(max(wtw,0._r8)+1.e-4_r8), 0._r8 ) 
             umf(k) = umf(k) * exp( 0.637_r8*(dpe/rho0j/g) * autodet )   
          end if      
          if(umf(k).eq.0._r8) wtw = -1._r8

          ! -------------------------------------- !
          ! Below block is just a dignostic output !
          ! -------------------------------------- ! 

          excessu_arr(k) = excessu
          excess0_arr(k) = excess0
          xc_arr(k) = xc
          aquad_arr(k) = aquad
          bquad_arr(k) = bquad
          cquad_arr(K) = cquad
          bogbot_arr(k) = bogbot
          bogtop_arr(k) = bogtop

          ! ------------------------------------------------------------------- !
          ! 'kbup' is the upper most layer in which cloud buoyancy  is positive ! 
          ! both at the base and top interface.  'kpen' is the upper most layer !
          ! up to cumulus can reach. Usually, 'kpen' is located higher than the !
          ! 'kbup'. Note we initialized these by 'kbup = krel' & 'kpen = krel'. !
          ! As explained before, it is possible that only 'kpen' is updated,    !
          ! while 'kbup' keeps its initialization value. For this case, current !
          ! scheme will simply turns-off penetrative entrainment fluxes and use ! 
          ! normal buoyancy-sorting fluxes for 'kbup <= k <= kpen-1' interfaces,!
          ! in order to describe shallow continental cumulus convection.        !
          ! ------------------------------------------------------------------- !
          
          if( bogbot .gt. 0._r8 .and. bogtop .gt. 0._r8 ) then 
             kbup = k
          end if
              
          if(wtw .le. 0._r8) then
             kpen = k
             go to 45
          end if

          wu(k) = sqrt(wtw)
          if(wu(k) .gt. 100._r8) then
           ! write(iulog,*) 'Too strong updraft'
             exit_wu(i) = 1._r8
             id_exit = .true.
             go to 333
          endif

          ! ---------------------------------------------------------------------------- !
          ! Iteration end due to 'rmaxfrac' constraint [ ***************************** ] ! 
          ! ---------------------------------------------------------------------------- !

          ! ---------------------------------------------------------------------- !
          ! Calculate updraft fractional area at the upper interface and set upper ! 
          ! limit to 'ufrc' by 'rmaxfrac'. In order to keep the consistency  among !
          ! ['ufrc','umf','wu (or wtw)'], if ufrc is limited by 'rmaxfrac', either !
          ! 'umf' or 'wu' should be changed. Although both 'umf' and 'wu (wtw)' at !
          ! the current upper interface are used for updating 'umf' & 'wu'  at the !
          ! next upper interface, 'umf' is a passive variable not influencing  the !
          ! buoyancy sorting process in contrast to 'wtw'. This is a reason why we !
          ! adjusted 'umf' instead of 'wtw'. In turn we updated 'fdr' here instead !
          ! of 'fer',  which guarantees  that all previously updated thermodynamic !
          ! variables at the upper interface before applying 'rmaxfrac' constraint !
          ! are already internally consistent,  even though 'ufrc'  is  limited by !
          ! 'rmaxfrac'. Thus, we don't need to go through interation loop again.If !
          ! If we update 'fer' however, we should go through above iteration loop. !
          ! ---------------------------------------------------------------------- !
            
          rhos0j = ps0(k) / ( r * 0.5_r8 * ( thv0bot(k+1) + thv0top(k) ) * exns0(k) )
          ufrc(k) = umf(k) / ( rhos0j * wu(k) )
          if(ufrc(k) .gt. rmaxfrac) then
             limit_ufrc(i) = 1._r8 
             ufrc(k) = rmaxfrac
             umf(k) = rmaxfrac * rhos0j * wu(k)
             fdr(k) = fer(k) - log( umf(k) / umf(km1) ) / dpe
          endif

          ! ------------------------------------------------------------ !
          ! Update environmental properties for at the mid-point of next !
          ! upper layer for use in buoyancy sorting.                     !
          ! ------------------------------------------------------------ ! 

          pe      = p0(k+1)
          dpe     = dp0(k+1)
          exne    = exn0(k+1)
          thvebot = thv0bot(k+1)
          thle    = thl0(k+1)
          qte     = qt0(k+1)
          ue      = u0(k+1)
          ve      = v0(k+1) 

       end do   ! End of cumulus updraft loop from the 'krel' layer to 'kpen' layer.
       
       ! ------------------------------------------------------------------------------- !
       ! Up to this point, we finished all of buoyancy sorting processes from the 'krel' !
       ! layer to 'kpen' layer: at the top interface of individual layers, we calculated !
       ! updraft and penetrative mass fluxes [ umf(k) & emf(k) = 0 ], updraft fractional !
       ! area [ ufrc(k) ],  updraft vertical velocity [ wu(k) ],  updraft  thermodynamic !
       ! variables [thlu(k),qtu(k),uu(k),vu(k),thvu(k)]. In the layer,we also calculated !
       ! fractional entrainment-detrainment rate [ fer(k), fdr(k) ], and detrainment ten !
       ! dency of water and ice from cumulus updraft [ dwten(k), diten(k) ]. In addition,!
       ! we updated and identified 'krel' and 'kpen' layer index, if any.  In the 'kpen' !
       ! layer, we calculated everything mentioned above except the 'wu(k)' and 'ufrc(k)'!
       ! since a real value of updraft vertical velocity is not defined at the kpen  top !
       ! interface (note 'ufrc' at the top interface of layer is calculated from 'umf(k)'!
       ! and 'wu(k)'). As mentioned before, special treatment is required when 'kbup' is !
       ! not updated and so 'kbup = krel'.                                               !
       ! ------------------------------------------------------------------------------- !
       
       ! ------------------------------------------------------------------------------ !
       ! During the 'iter_scaleh' iteration loop, non-physical ( with non-zero values ) !
       ! values can remain in the variable arrays above (also 'including' in case of wu !
       ! and ufrc at the top interface) the 'kpen' layer. This can happen when the kpen !
       ! layer index identified from the 'iter_scaleh = 1' iteration loop is located at !
       ! above the kpen layer index identified from   'iter_scaleh = 3' iteration loop. !
       ! Thus, in the following calculations, we should only use the values in each     !
       ! variables only up to finally identified 'kpen' layer & 'kpen' interface except ! 
       ! 'wu' and 'ufrc' at the top interface of 'kpen' layer.    Note that in order to !
       ! prevent any problems due to these non-physical values, I re-initialized    the !
       ! values of [ umf(kpen:mkx), emf(kpen:mkx), dwten(kpen+1:mkx), diten(kpen+1:mkx),! 
       ! fer(kpen:mkx), fdr(kpen+1:mkx), ufrc(kpen:mkx) ] to be zero after 'iter_scaleh'!
       ! do loop.                                                                       !
       ! ------------------------------------------------------------------------------ !
       
 45    continue

       ! ----------------------------- !
       ! Set fer(kpen) to zero or not. !
       ! ----------------------------- !

       if( set_zeroferkpen ) then      

          thlu(kpen) = thlu(kpen-1)
          qtu(kpen) = qtu(kpen-1)       
          call conden(ps0(kpen),thlu(kpen),qtu(kpen),thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
          end if
          if( (qlj + qij) .gt. criqc ) then
            exql = ( ( qlj + qij ) - criqc ) * qlj / ( qlj + qij )
            exqi = ( ( qlj + qij ) - criqc ) * qij / ( qlj + qij )
            qtu(kpen)  = qtu(kpen) - exql - exqi
            thlu(kpen) = thlu(kpen) + (xlv/cp/exns0(kpen))*exql + (xls/cp/exns0(kpen))*exqi 
          endif      
          call conden(ps0(kpen),thlu(kpen),qtu(kpen),thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
          end if
          thvu(kpen) = thj * ( 1._r8 + zvir * qvj - qlj - qij )
          bogbot = rbuoy * ( thvu(kpen-1)/thv0bot(kpen) - 1._r8 )
          bogtop = rbuoy * ( thvu(kpen)/thv0top(kpen)   - 1._r8 )
          delbog = bogtop - bogbot
          drage = 0._r8
 
          fer(kpen) = 0._r8 
          fdr(kpen) = rei(kpen)
          xc_arr(kpen) = 0._r8
          bogbot_arr(kpen) = bogbot
          bogtop_arr(kpen) = bogtop

       endif

       ! ------------------------------------------------------------------------------ !
       ! Calculate 'ppen( < 0 )', updarft penetrative distance from the lower interface !
       ! of 'kpen' layer. Note that bogbot & bogtop at the 'kpen' layer either when fer !
       ! is zero or non-zero was already calculated above.                              !
       ! It seems that below qudarature solving formula is valid only when bogbot < 0.  !
       ! Below solving equation is clearly wrong ! I should revise this !               !
       ! ------------------------------------------------------------------------------ ! 
            
       if(drage.eq.0._r8) then
          aquad =  ( bogtop - bogbot ) / ( ps0(kpen) - ps0(kpen-1) )
          bquad =  2._r8 * bogbot
          cquad = -wu(kpen-1)**2 * rho0j
          call roots(aquad,bquad,cquad,xc1,xc2,status)
          if( status .eq. 0 ) then
              if( xc1 .le. 0._r8 .and. xc2 .le. 0._r8 ) then
                  ppen = max( xc1, xc2 )
                  ppen = min( 0._r8,max( -dp0(kpen), ppen ) )  
              elseif( xc1 .gt. 0._r8 .and. xc2 .gt. 0._r8 ) then
                  ppen = -dp0(kpen)
                  write(iulog,*) 'Warning : UW-Cumulus penetrates upto kpen interface'
              else
                  ppen = min( xc1, xc2 )
                  ppen = min( 0._r8,max( -dp0(kpen), ppen ) )  
              endif
          else
              ppen = -dp0(kpen)
              write(iulog,*) 'Warning : UW-Cumulus penetrates upto kpen interface'
          endif       
       else 
          ppen = compute_ppen(wtwb,drage,bogbot,bogtop,rho0j,dp0(kpen))
       endif
       if( ppen.eq.-dp0(kpen).or.ppen.eq.0._r8 ) limit_ppen(i) = 1._r8

       ! ----------------------------------------------------------------- !
       ! Re-calculate the amount of expelled condensate from cloud updraft !
       ! at the cumulus top. This is necessary for refined calculations of !
       ! bulk cloud microphysics at the cumulus top. Note that ppen < 0._r8   !
       ! In the below, I explicitly calculate 'thlu_top' & 'qtu_top' by    !
       ! using non-zero 'fer(kpen)'.                                       !    
       ! ----------------------------------------------------------------- !

       if( fer(kpen)*(-ppen) .lt. 1.e-4_r8 ) then
           thlu_top = thlu(kpen-1) + & 
                    ( thl0(kpen) + ssthl0(kpen) * (-ppen) / 2._r8 - thlu(kpen-1) ) * fer(kpen) * (-ppen)
           qtu_top  =  qtu(kpen-1) + & 
                     ( qt0(kpen) + ssqt0(kpen) * (-ppen) / 2._r8  -  qtu(kpen-1) ) * fer(kpen) * (-ppen)
       else
           thlu_top = ( thl0(kpen) + ssthl0(kpen) / fer(kpen) - ssthl0(kpen) * (-ppen) / 2._r8 ) - &
                      ( thl0(kpen) + ssthl0(kpen) * (-ppen) / 2._r8 - thlu(kpen-1) + ssthl0(kpen) / fer(kpen) ) * & 
                        exp(-fer(kpen) * (-ppen))
           qtu_top  = ( qt0(kpen)  + ssqt0(kpen) / fer(kpen) - ssqt0(kpen) * (-ppen) / 2._r8 ) - &  
                      ( qt0(kpen)  + ssqt0(kpen) * (-ppen) / 2._r8 - qtu(kpen-1) + ssqt0(kpen) / fer(kpen) ) * &
                        exp(-fer(kpen) * (-ppen))
       end if

       call conden(ps0(kpen-1)+ppen,thlu_top,qtu_top,thj,qvj,qlj,qij,qse,id_check,qsat)
       if(id_check.eq.1) then
          exit_conden(i) = 1._r8
          id_exit = .true.
          go to 333
       end if
       exntop = ((ps0(kpen-1)+ppen)/p00)**rovcp
       if( (qlj + qij) .gt. criqc ) then
            dwten(kpen) = ( ( qlj + qij ) - criqc ) * qlj / ( qlj + qij )
            diten(kpen) = ( ( qlj + qij ) - criqc ) * qij / ( qlj + qij )
            qtu_top  = qtu_top - dwten(kpen) - diten(kpen)
            thlu_top = thlu_top + (xlv/cp/exntop)*dwten(kpen) + (xls/cp/exntop)*diten(kpen) 
       else
            dwten(kpen) = 0._r8
            diten(kpen) = 0._r8
       endif      
 
       ! ----------------------------------------------------------------------- !
       ! Calculate cumulus scale height as the top height that cumulus can reach.!
       ! ----------------------------------------------------------------------- !
       
       rhos0j = ps0(kpen-1)/(r*0.5_r8*(thv0bot(kpen)+thv0top(kpen-1))*exns0(kpen-1))  
       cush   = zs0(kpen-1) - ppen/rhos0j/g
       scaleh = cush 

    end do   ! End of 'iter_scaleh' loop.   


       ! -------------------------------------------------------------------- !   
       ! The 'shallowCu' is logical identifier saying whether cumulus updraft !
       ! overcome the buoyancy barrier just above the PBL top. If it is true, !
       ! cumulus did not overcome the barrier -  this is a shallow convection !
       ! with negative cloud buoyancy, mimicking  shallow continental cumulus !
       ! convection. Depending on 'shallowCu' parameter, treatment of heat  & !
       ! moisture fluxes at the entraining interfaces, 'kbup <= k < kpen - 1' !
       ! will be set up in a different ways, as will be shown later.          !
       ! -------------------------------------------------------------------- !
 
       if( kbup.eq.krel ) then 
           shallowCu =.true.
           limit_shcu(i) = 1._r8
       else
           shallowCu =.false.
           limit_shcu(i) = 0._r8
       endif  
       
       ! ------------------------------------------------------------------ !
       ! Filtering of unerasonable cumulus adjustment here.  This is a very !
       ! important process which should be done cautiously. Various ways of !
       ! filtering are possible depending on cases mainly using the indices !
       ! of key layers - 'klcl','kinv','krel','klfc','kbup','kpen'. At this !
       ! stage, the followings are all possible : 'kinv >= 2', 'klcl >= 1', !
       ! 'krel >= kinv', 'kbup >= krel', 'kpen >= krel'. I must design this !
       ! filtering very cautiously, in such that none of  realistic cumulus !
       ! convection is arbitrarily turned-off. Potentially, I might turn-off! 
       ! cumulus convection if layer-mean 'ql > 0' in the 'kinv-1' layer,in !
       ! order to suppress cumulus convection growing, based at the Sc top. ! 
       ! This is one of potential future modifications. Note that ppen < 0. !
       ! ------------------------------------------------------------------ !

       cldhgt = ps0(kpen-1) + ppen
       if( shallowCu ) then
           ! write(iulog,*) 'shallowCu - did not overcome initial buoyancy barrier'
           exit_cufilter(i) = 1._r8
           id_exit = .true.
           go to 333
       end if
       ! Limit 'additional shallow cumulus' for DYCOMS simulation.
       ! if( cldhgt.ge.88000._r8 ) then
       !     id_exit = .true.
       !     go to 333
       ! end if
       
       ! ------------------------------------------------------------------------------ !
       ! Re-initializing some key variables above the 'kpen' layer in order to suppress !
       ! the influence of non-physical values above 'kpen', in association with the use !
       ! of 'iter_scaleh' loop. Note that umf, emf,  ufrc are defined at the interfaces !
       ! (0:mkx), while 'dwten','diten', 'fer', 'fdr' are defined at layer mid-points.  !
       ! Initialization of 'fer' and 'fdr' is for correct writing purpose of diagnostic !
       ! output. Note that we set umf(kpen)=emf(kpen)=ufrc(kpen)=0, in consistent  with !
       ! wtw < 0  at the top interface of 'kpen' layer. However, we still have non-zero !
       ! expelled cloud condensate in the 'kpen' layer.                                 !
       ! ------------------------------------------------------------------------------ !

       umf(kpen:mkx) = 0._r8
       emf(kpen:mkx) = 0._r8
       ufrc(kpen:mkx) = 0._r8
       dwten(kpen+1:mkx) = 0._r8
       diten(kpen+1:mkx) = 0._r8
       fer(kpen+1:mkx) = 0._r8
       fdr(kpen+1:mkx) = 0._r8
       ! fer(kpen) = 0._r8 
       ! fdr(kpen) = rei(kpen)
       
       ! ------------------------------------------------------------------------ !
       ! Calculate downward penetrative entrainment mass flux, 'emf(k) < 0',  and !
       ! thermodynamic properties of penetratively entrained airs at   entraining !
       ! interfaces. emf(k) is defined from the top interface of the  layer  kbup !
       ! to the bottom interface of the layer 'kpen'. Note even when  kbup = krel,!
       ! i.e.,even when 'kbup' was not updated in the above buoyancy  sorting  do !
       ! loop (i.e., 'kbup' remains as the initialization value),   below do loop !
       ! of penetrative entrainment flux can be performed without  any conceptual !
       ! or logical problems, because we have already computed all  the variables !
       ! necessary for performing below penetrative entrainment block.            !
       ! In the below 'do' loop, 'k' is an interface index at which non-zero 'emf'! 
       ! (penetrative entrainment mass flux) is calculated. Since cumulus updraft !
       ! is negatively buoyant in the layers between the top interface of 'kbup'  !
       ! layer (interface index, kbup) and the top interface of 'kpen' layer, the !
       ! fractional lateral entrainment, fer(k) within these layers will be close !
       ! to zero - so it is likely that only strong lateral detrainment occurs in !
       ! thses layers. Under this situation,we can easily calculate the amount of !
       ! detrainment cumulus air into these negatively buoyanct layers by  simply !
       ! comparing cumulus updraft mass fluxes between the base and top interface !
       ! of each layer: emf(k) = emf(k-1)*exp(-fdr(k)*dp0(k))                     !
       !                       ~ emf(k-1)*(1-rei(k)*dp0(k))                       !
       !                emf(k-1)-emf(k) ~ emf(k-1)*rei(k)*dp0(k)                  !
       ! Current code assumes that about 'rpen~10' times of these detrained  mass !
       ! are penetratively re-entrained down into the 'k-1' interface. And all of !
       ! these detrained masses are finally dumped down into the top interface of !
       ! 'kbup' layer. Thus, the amount of penetratively entrained air across the !
       ! top interface of 'kbup' layer with 'rpen~10' becomes too large.          !
       ! Note that this penetrative entrainment part can be completely turned-off !
       ! and we can simply use normal buoyancy-sorting involved turbulent  fluxes !
       ! by modifying 'penetrative entrainment fluxes' part below.                !
       ! ------------------------------------------------------------------------ !
       
       ! -----------------------------------------------------------------------!
       ! Calculate entrainment mass flux and conservative scalars of entraining !
       ! free air at interfaces of 'kbup <= k < kpen - 1'                       !
       ! ---------------------------------------------------------------------- !
 
       do k = 0, mkx
          thlu_emf(k) = thlu(k)
          qtu_emf(k)  = qtu(k)
          uu_emf(k)   = uu(k)
          vu_emf(k)   = vu(k)
       end do

       do k = kpen - 1, kbup, -1  ! Here, 'k' is an interface index at which
                                  ! penetrative entrainment fluxes are calculated. 
                                  
          rhos0j = ps0(k) / ( r * 0.5_r8 * ( thv0bot(k+1) + thv0top(k) ) * exns0(k) )

          if( k .eq. kpen - 1 ) then

             ! ------------------------------------------------------------------------ ! 
             ! Note that 'ppen' has already been calculated in the above 'iter_scaleh'  !
             ! loop assuming zero lateral entrainmentin the layer 'kpen'.               !
             ! ------------------------------------------------------------------------ !       
             
             ! -------------------------------------------------------------------- !
             ! Calculate returning mass flux, emf ( < 0 )                           !
             ! Current penetrative entrainment rate with 'rpen~10' is too large and !
             ! future refinement is necessary including the definition of 'thl','qt'! 
             ! of penetratively entrained air.  Penetratively entrained airs across !
             ! the 'kpen-1' interface is assumed to have the properties of the base !
             ! interface of 'kpen' layer. Note that 'emf ~ - umf/ufrc = - w * rho'. !
             ! Thus, below limit sets an upper limit of |emf| to be ~ 10cm/s, which !
             ! is very loose constraint. Here, I used more restricted constraint on !
             ! the limit of emf, assuming 'emf' cannot exceed a net mass within the !
             ! layer above the interface. Similar to the case of warming and drying !
             ! due to cumulus updraft induced compensating subsidence,  penetrative !
             ! entrainment induces compensating upwelling -     in order to prevent !  
             ! numerical instability in association with compensating upwelling, we !
             ! should similarily limit the amount of penetrative entrainment at the !
             ! interface by the amount of masses within the layer just above the    !
             ! penetratively entraining interface.                                  !
             ! -------------------------------------------------------------------- !
             
             if((umf(k)*ppen*rei(kpen)*rpen).lt.-0.1*rhos0j)         limit_emf(i) = 1._r8
             if((umf(k)*ppen*rei(kpen)*rpen).lt.-0.9*dp0(kpen)/g/dt) limit_emf(i) = 1._r8             

             emf(k)  = max(max(umf(k)*ppen*rei(kpen)*rpen, -0.1_r8*rhos0j), -0.9_r8*dp0(kpen)/g/dt)
             thlu_emf(k) = thl0(kpen) + ssthl0(kpen) * ( ps0(k) - p0(kpen) )
             qtu_emf(k)  = qt0(kpen)  + ssqt0(kpen)  * ( ps0(k) - p0(kpen) )
             uu_emf(k)   = u0(kpen)   + ssu0(kpen)   * ( ps0(k) - p0(kpen) )     
             vu_emf(k)   = v0(kpen)   + ssv0(kpen)   * ( ps0(k) - p0(kpen) )   

          else ! if(k.lt.kpen-1). 
              
             ! --------------------------------------------------------------------------- !
             ! Note we are coming down from the higher interfaces to the lower interfaces. !
             ! Also note that 'emf < 0'. So, below operation is a summing not subtracting. !
             ! In order to ensure numerical stability, I imposed a modified correct limit  ! 
             ! of '-0.9*dp0(k+1)/g/dt' on emf(k).                                          !
             ! --------------------------------------------------------------------------- !

             if((emf(k+1)-umf(k)*dp0(k+1)*rei(k+1)*rpen).lt.-0.1_r8*rhos0j)        limit_emf(i) = 1
             if((emf(k+1)-umf(k)*dp0(k+1)*rei(k+1)*rpen).lt.-0.9_r8*dp0(k+1)/g/dt) limit_emf(i) = 1           ! new code
           ! if((emf(k+1)-umf(k)*dp0(k+1)*rei(k+1)*rpen).lt.(emf(k+1)-0.9_r8*dp0(k+1)/g/dt)) limit_emf(i) = 1 ! old code.
         
           ! emf(k)  = max(max(emf(k+1)-umf(k)*dp0(k+1)*rei(k+1)*rpen, -0.1_r8*rhos0j), emf(k+1)-0.9_r8*dp0(k+1)/g/dt) ! old code.
             emf(k)  = max(max(emf(k+1)-umf(k)*dp0(k+1)*rei(k+1)*rpen, -0.1_r8*rhos0j), -0.9_r8*dp0(k+1)/g/dt )        ! new code
             if( abs(emf(k)) .gt. abs(emf(k+1)) ) then
                thlu_emf(k) = ( thlu_emf(k+1) * emf(k+1) + thl0(k+1) * ( emf(k) - emf(k+1) ) ) / emf(k)
                qtu_emf(k)  = ( qtu_emf(k+1)  * emf(k+1) + qt0(k+1)  * ( emf(k) - emf(k+1) ) ) / emf(k)
                uu_emf(k)   = ( uu_emf(k+1)   * emf(k+1) + u0(k+1)   * ( emf(k) - emf(k+1) ) ) / emf(k)
                vu_emf(k)   = ( vu_emf(k+1)   * emf(k+1) + v0(k+1)   * ( emf(k) - emf(k+1) ) ) / emf(k)
             else   
                thlu_emf(k) = thl0(k+1)
                qtu_emf(k)  =  qt0(k+1)
                uu_emf(k)   =   u0(k+1)
                vu_emf(k)   =   v0(k+1)
             endif   
                     
          endif

          ! ---------------------------------------------------------------------------- !
          ! In this GCM modeling framework,  all what we should do is to calculate  heat !
          ! and moisture fluxes at the given geometrically-fixed height interfaces -  we !
          ! don't need to worry about movement of material height surface in association !
          ! with compensating subsidence or unwelling, in contrast to the bulk modeling. !
          ! In this geometrically fixed height coordinate system, heat and moisture flux !
          ! at the geometrically fixed height handle everything - a movement of material !
          ! surface is implicitly treated automatically. Note that in terms of turbulent !
          ! heat and moisture fluxes at model interfaces, both the cumulus updraft  mass !
          ! flux and penetratively entraining mass flux play the same role -both of them ! 
          ! warms and dries the 'kbup' layer, cools and moistens the 'kpen' layer,   and !
          ! cools and moistens any intervening layers between 'kbup' and 'kpen' layers.  !
          ! It is important to note these identical roles on turbulent heat and moisture !
          ! fluxes of 'umf' and 'emf'.                                                   !
          ! When 'kbup' is a stratocumulus-topped PBL top interface,  increase of 'rpen' !
          ! is likely to strongly diffuse stratocumulus top interface,  resulting in the !
          ! reduction of cloud fraction. In this sense, the 'kbup' interface has a  very !
          ! important meaning and role : across the 'kbup' interface, strong penetrative !
          ! entrainment occurs, thus any sharp gradient properties across that interface !
          ! are easily diffused through strong mass exchange. Thus, an initialization of ! 
          ! 'kbup' (and also 'kpen') should be done very cautiously as mentioned before. ! 
          ! In order to prevent this stron diffusion for the shallow cumulus convection  !
          ! based at the Sc top, it seems to be good to initialize 'kbup = krel', rather !
          ! that 'kbup = krel-1'.                                                        !
          ! ---------------------------------------------------------------------------- !
          
       end do

       !------------------------------------------------------------------ !
       !                                                                   ! 
       ! Compute turbulent heat, moisture, momentum flux at all interfaces !
       !                                                                   !
       !------------------------------------------------------------------ !
       ! It is very important to note that in calculating turbulent fluxes !
       ! below, we must not double count turbulent flux at any interefaces.!
       ! In the below, turbulent fluxes at the interfaces (interface index !
       ! k) are calculated by the following 4 blocks in consecutive order: !
       !                                                                   !
       ! (1) " 0 <= k <= kinv - 1 "  : PBL fluxes.                         !
       !     From 'fluxbelowinv' using reconstructed PBL height. Currently,!
       !     the reconstructed PBLs are independently calculated for  each !
       !     individual conservative scalar variables ( qt, thl, u, v ) in !
       !     each 'fluxbelowinv',  instead of being uniquely calculated by !
       !     using thvl. Turbulent flux at the surface is assumed to be 0. !
       ! (2) " kinv <= k <= krel - 1 " : Non-buoyancy sorting fluxes       !
       !     Assuming cumulus mass flux  and cumulus updraft thermodynamic !
       !     properties (except u, v which are modified by the PGFc during !
       !     upward motion) are conserved during a updraft motion from the !
       !     PBL top interface to the release level. If these layers don't !
       !     exist (e,g, when 'krel = kinv'), then  current routine do not !
       !     perform this routine automatically. So I don't need to modify !
       !     anything.                                                     ! 
       ! (3) " krel <= k <= kbup - 1 " : Buoyancy sorting fluxes           !
       !     From laterally entraining-detraining buoyancy sorting plumes. ! 
       ! (4) " kbup <= k < kpen-1 " : Penetrative entrainment fluxes       !
       !     From penetratively entraining plumes,                         !
       !                                                                   !
       ! In case of normal situation, turbulent interfaces  in each groups !
       ! are mutually independent of each other. Thus double flux counting !
       ! or ambiguous flux counting requiring the choice among the above 4 !
       ! groups do not occur normally. However, in case that cumulus plume !
       ! could not completely overcome the buoyancy barrier just above the !
       ! PBL top interface and so 'kbup = krel' (.shallowCu=.true.) ( here,!
       ! it can be either 'kpen = krel' as the initialization, or ' kpen > !
       ! krel' if cumulus updraft just penetrated over the top of  release !
       ! layer ). If this happens, we should be very careful in organizing !
       ! the sequence of the 4 calculation routines above -  note that the !
       ! routine located at the later has the higher priority.  Additional ! 
       ! feature I must consider is that when 'kbup = kinv - 1' (this is a !
       ! combined situation of 'kbup=krel-1' & 'krel = kinv' when I  chose !
       ! 'kbup=krel-1' instead of current choice of 'kbup=krel'), a strong !
       ! penetrative entrainment fluxes exists at the PBL top interface, & !
       ! all of these fluxes are concentrated (deposited) within the layer ! 
       ! just below PBL top interface (i.e., 'kinv-1' layer). On the other !
       ! hand, in case of 'fluxbelowinv', only the compensating subsidence !
       ! effect is concentrated in the 'kinv-1' layer and 'pure' turbulent !
       ! heat and moisture fluxes ( 'pure' means the fluxes not associated !
       ! with compensating subsidence) are linearly distributed throughout !
       ! the whole PBL. Thus different choice of the above flux groups can !
       ! produce very different results. Output variable should be written !
       ! consistently to the choice of computation sequences.              !
       ! When the case of 'kbup = krel(-1)' happens,another way to dealing !
       ! with this case is to simply ' exit ' the whole cumulus convection !
       ! calculation without performing any cumulus convection.     We can !
       ! choose this approach by specifying a condition in the  'Filtering !
       ! of unreasonable cumulus adjustment' just after 'iter_scaleh'. But !
       ! this seems not to be a good choice (although this choice was used !
       ! previous code ), since it might arbitrary damped-out  the shallow !
       ! cumulus convection over the continent land, where shallow cumulus ! 
       ! convection tends to be negatively buoyant.                        !
       ! ----------------------------------------------------------------- !  

       ! --------------------------------------------------- !
       ! 1. PBL fluxes :  0 <= k <= kinv - 1                 !
       !    All the information necessary to reconstruct PBL ! 
       !    height are passed to 'fluxbelowinv'.             !
       ! --------------------------------------------------- !

       xsrc  = qtsrc
       xmean = qt0(kinv)
       xtop  = qt0(kinv+1) + ssqt0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
       xbot  = qt0(kinv-1) + ssqt0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )        
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
       qtflx(0:kinv-1) = xflx(0:kinv-1)

       xsrc  = thlsrc
       xmean = thl0(kinv)
       xtop  = thl0(kinv+1) + ssthl0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
       xbot  = thl0(kinv-1) + ssthl0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )        
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
       slflx(0:kinv-1) = cp * exns0(0:kinv-1) * xflx(0:kinv-1)

       xsrc  = usrc
       xmean = u0(kinv)
       xtop  = u0(kinv+1) + ssu0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
       xbot  = u0(kinv-1) + ssu0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
       uflx(0:kinv-1) = xflx(0:kinv-1)

       xsrc  = vsrc
       xmean = v0(kinv)
       xtop  = v0(kinv+1) + ssv0(kinv+1) * ( ps0(kinv)   - p0(kinv+1) )
       xbot  = v0(kinv-1) + ssv0(kinv-1) * ( ps0(kinv-1) - p0(kinv-1) )
       call fluxbelowinv( cbmf, ps0(0:mkx), mkx, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
       vflx(0:kinv-1) = xflx(0:kinv-1)

       ! -------------------------------------------------------------- !
       ! 2. Non-buoyancy sorting fluxes : kinv <= k <= krel - 1         !
       !    Note that when 'krel = kinv', below block is never executed !
       !    as in a desirable, expected way ( but I must check  if this !
       !    is the case ). The non-buoyancy sorting fluxes are computed !
       !    only when 'krel > kinv'.                                    !
       ! -------------------------------------------------------------- !          

       uplus = 0._r8
       vplus = 0._r8
       do k = kinv, krel - 1 
          qtflx(k) = cbmf * ( qtsrc  - (  qt0(k+1) +  ssqt0(k+1) * ( ps0(k) - p0(k+1) ) ) )          
          slflx(k) = cbmf * ( thlsrc - ( thl0(k+1) + ssthl0(k+1) * ( ps0(k) - p0(k+1) ) ) ) * cp * exns0(k)
          uplus = uplus + PGFc * ssu0(k) * ( ps0(k) - ps0(k-1) )
          vplus = vplus + PGFc * ssv0(k) * ( ps0(k) - ps0(k-1) )
          uflx(k) = cbmf * ( usrc + uplus -  (  u0(k+1)  +   ssu0(k+1) * ( ps0(k) - p0(k+1) ) ) ) 
          vflx(k) = cbmf * ( vsrc + vplus -  (  v0(k+1)  +   ssv0(k+1) * ( ps0(k) - p0(k+1) ) ) )
       end do

       ! ------------------------------------------------------------------------ !
       ! 3. Buoyancy sorting fluxes : krel <= k <= kbup - 1                       !
       !    In case that 'kbup = krel - 1 ' ( or even in case 'kbup = krel' ),    ! 
       !    buoyancy sorting fluxes are not calculated, which is consistent,      !
       !    desirable feature.                                                    !  
       ! ------------------------------------------------------------------------ !

       do k = krel, kbup - 1      
          kp1 = k + 1
          slflx(k) = cp * exns0(k) * umf(k)*(thlu(k) - (thl0(kp1) + ssthl0(kp1)*(ps0(k) - p0(kp1))))
          qtflx(k) = umf(k)*(qtu(k) - (qt0(kp1) + ssqt0(kp1)*(ps0(k) - p0(kp1))))
          uflx(k) = umf(k)*(uu(k) - (u0(kp1) + ssu0(kp1)*(ps0(k) - p0(kp1))))
          vflx(k) = umf(k)*(vu(k) - (v0(kp1) + ssv0(kp1)*(ps0(k) - p0(kp1))))
       end do

       ! ------------------------------------------------------------------------- !
       ! 4. Penetrative entrainment fluxes : kbup <= k <= kpen - 1                 !
       !    The only confliction that can happen is when 'kbup = kinv-1'. For this !
       !    case, turbulent flux at kinv-1 is calculated  both from 'fluxbelowinv' !
       !    and here as penetrative entrainment fluxes.  Since penetrative flux is !
       !    calculated later, flux at 'kinv - 1 ' will be that of penetrative flux.!
       !    However, turbulent flux calculated at 'kinv - 1' from penetrative entr.!
       !    is less attractable,  since more reasonable turbulent flux at 'kinv-1' !
       !    should be obtained from 'fluxbelowinv', by considering  re-constructed ! 
       !    inversion base height. This conflicting problem can be solved if we can!
       !    initialize 'kbup = krel', instead of kbup = krel - 1. This choice seems!
       !    to be more reasonable since it is not conflicted with 'fluxbelowinv' in!
       !    calculating fluxes at 'kinv - 1' ( for this case, flux at 'kinv-1' is  !
       !    always from 'fluxbelowinv' ), and flux at 'krel-1' is calculated from  !
       !    the non-buoyancy sorting flux without being competed with penetrative  !
       !    entrainment fluxes. Even when we use normal cumulus flux instead of    !
       !    penetrative entrainment fluxes at 'kbup <= k <= kpen-1' interfaces,    !
       !    the initialization of kbup=krel perfectly works without any conceptual !
       !    confliction. Thus it seems to be much better to choose 'kbup = krel'   !
       !    initialization of 'kbup', which is current choice.                     !
       !    Note that below formula uses conventional updraft cumulus fluxes for   !
       !    shallow cumulus which did not overcome the first buoyancy barrier above!
       !    PBL top while uses penetrative entrainment fluxes for the other cases  !
       !    'kbup <= k <= kpen-1' interfaces. Depending on cases, however, I can   !
       !    selelct different choice.                                              !
       ! --------------------------------------------------------------------------------------------------- !
       !  slflx(k) = cp * exns0(k) * umf(k)*(thlu(k)     - (thl0(kp1) + ssthl0(kp1)*(ps0(k) - p0(kp1)))) + & !
       !             cp * exns0(k) * emf(k)*(thlu_emf(k) - (thl0(k)   + ssthl0(k)  *(ps0(k) - p0(k))))       !
       !  qtflx(k) = umf(k)*(qtu(k)     - (qt0(kp1) + ssqt0(kp1)*(ps0(k) - p0(kp1)))) +                    & !
       !             emf(k)*(qtu_emf(k) - (qt0(k)   + ssqt0(k)  *(ps0(k) - p0(k))))                          !
       !  uflx(k) = umf(k)*(uu(k)     - (u0(kp1) + ssu0(kp1)*(ps0(k) - p0(kp1)))) +                        & !
       !             emf(k)*(uu_emf(k) - (u0(k)   + ssu0(k)  *(ps0(k) - p0(k))))                             ! 
       !  vflx(k) = umf(k)*(vu(k)     - (v0(kp1) + ssv0(kp1)*(ps0(k) - p0(kp1)))) +                        & !
       !             emf(k)*(vu_emf(k) - (v0(k)   + ssv0(k)  *(ps0(k) - p0(k))))                             !
       ! --------------------------------------------------------------------------------------------------- !

       do k = kbup, kpen - 1      
          kp1 = k + 1
          if( shallowCu ) then
              slflx(k) = cp * exns0(k) * umf(k)*(thlu(k) - (thl0(kp1) + ssthl0(kp1)*(ps0(k) - p0(kp1))))
              qtflx(k) = umf(k)*(qtu(k) - (qt0(kp1) + ssqt0(kp1)*(ps0(k) - p0(kp1))))
              uflx(k) = umf(k)*(uu(k)  - (u0(kp1) + ssu0(kp1)*(ps0(k) - p0(kp1))))
              vflx(k) = umf(k)*(vu(k)  - (v0(kp1) + ssv0(kp1)*(ps0(k) - p0(kp1))))
          else
              slflx(k) = cp * exns0(k) * emf(k)*(thlu_emf(k) - (thl0(k) + ssthl0(k)  *(ps0(k) - p0(k))))
              qtflx(k) = emf(k)*(qtu_emf(k) - (qt0(k)  + ssqt0(k) *(ps0(k) - p0(k)))) 
              uflx(k) = emf(k)*(uu_emf(k)  - (u0(k)   + ssu0(k)  *(ps0(k) - p0(k)))) 
              vflx(k) = emf(k)*(vu_emf(k)  - (v0(k)   + ssv0(k)  *(ps0(k) - p0(k)))) 
          endif 
       end do

       ! ------------------------------------------- !
       ! Turn-off cumulus momentum flux as an option !
       ! ------------------------------------------- !

       if(.not.use_momenflx) then
         uflx(0:mkx) = 0._r8
         vflx(0:mkx) = 0._r8
       endif       

       ! --------------------------------------------------------------------------------- !
       ! Calculate net upward mass flux (uemf(k)) by summing "umf(k)"(>0) and "emf(k)"(<0) !
       ! at the interfaces "kinv-1 <= k <= kpen-1".  Set to zero for the other interfaces. !   
       ! Also, calculate net downward compensating subsidence ( comsub > 0 if compensating ! 
       ! subsidence is downward, or if uemf is positive ) at  a layer mid-point for layers !
       ! 'kinv <= k <= kpen". Set to zero for the other layers. Note that 'uemf(kpen) = 0' !
       ! in calculating compensating using the below formula is completely OK.             !
       ! --------------------------------------------------------------------------------- !
       
       uemf(0:mkx) = 0._r8
       uemf(kinv-1:krel-1) = cbmf
       uemf(krel:kbup-1) = umf(krel:kbup-1)
       uemf(kbup:kpen-1) = umf(kbup:kpen-1) + emf(kbup:kpen-1)

       comsub(1:mkx) = 0._r8
       do k = kinv, kpen
          comsub(k) = 0.5_r8 * ( uemf(k) + uemf(k-1) ) 
       end do    

       ! ---------------------------------------------------------------------------------- !
       ! Calculate condensate sink term induced by compensating subsidence   in association !  
       ! with cumulus updraft and penetrative entrainment term. Use centered difference for !
       ! layers " kinv + 1 <= k <= kpen - 1" while use upward difference for the  ambiguous ! 
       ! invesion layer, kinv. Note that "qcten_sink(k) > 0" means "qcten > 0", i.e.,       !
       ! condensate is generated by compensating subsidence. "cldfrct(k)" is total cloud    !
       ! fraction retrived from the physical buffer of previous time step.                  !
       ! Q : should I do this in the penetrative entrainment layer, kbup+1<=k<kpen-1(kpen)? !
       ! I am pretty sure that this is double counting since the effects of compensting sub !
       ! sidence or upwelling are already fully incorporated in the calculation of turbulent!
       ! heat ( slflx ) and moisture ( qtflx ) fluxes. Thus, we should not include this sink!
       ! tenddency term of condensate in the following equations.                           !
       ! ---------------------------------------------------------------------------------- !
       
       qcten_sink(:mkx) = 0._r8
       qlten_sink(:mkx) = 0._r8
       qiten_sink(:mkx) = 0._r8

       qcten_sink(kinv) = g * comsub(kinv) * ( cldfrct(kinv) - concldfrct(kinv) ) * & 
                             ( qs0(kinv+1) - qs0(kinv) ) / ( p0(kinv) - p0(kinv+1) )
       qcten_sink(kinv) = min( 0._r8, qcten_sink(kinv) )

       do k = kinv + 1, kpen - 1
         qcten_sink(k) = g * comsub(k) * ( cldfrct(k) - concldfrct(k) ) * &
                                         ( qs0(k+1) - qs0(k-1) ) / ( p0(k-1) - p0(k+1) )
         qcten_sink(k) = min( 0._r8, qcten_sink(k) )
       end do

       do k = kinv, kpen - 1  
         if( ( ql0(k) + qi0(k) ) .gt. 0._r8 ) then
           qlten_sink(k) = qcten_sink(k) * ql0(k) / ( ql0(k) + qi0(k) )
           qiten_sink(k) = qcten_sink(k) * qi0(k) / ( ql0(k) + qi0(k) )
         else
           qlten_sink(k) = 0._r8
           qiten_sink(k) = 0._r8
         end if 
       end do

       ! --------------------------------------------- !
       !                                               !
       ! Calculate convective tendencies at each layer ! 
       !                                               !
       ! --------------------------------------------- !
       
       ! ----------------- !
       ! Momentum tendency !
       ! ----------------- !
       
       do k = 1, kpen
          km1 = k - 1 
          uten(k) = ( uflx(km1) - uflx(k) ) * g / dp0(k)
          vten(k) = ( vflx(km1) - vflx(k) ) * g / dp0(k) 
          uf(k)   = u0(k) + uten(k) * dt
          vf(k)   = v0(k) + vten(k) * dt
       end do        

       ! ----------------------------------------------------------------- !
       ! Tendencies of thermodynamic variables.                            ! 
       ! This part requires a careful treatment of bulk cloud microphysics.!
       ! Relocations of 'precipitable condensates' either into the surface ! 
       ! or into the tendency of 'krel' layer will be performed just after !
       ! finishing the below 'do-loop'.                                    !        
       ! ----------------------------------------------------------------- !
       
       rliq    = 0._r8
       rainflx = 0._r8
       snowflx = 0._r8

       do k = 1, kpen

          km1 = k - 1

          ! ------------------------------------------------------------------------------ !
          ! Compute 'slten', 'qtten', 'qvten', 'qlten', 'qiten', and 'sten'                !
          !                                                                                !
          ! Key assumptions made in this 'cumulus scheme' are :                            !
          ! 1. Cumulus updraft expels condensate into the environment at the top interface !
          !    of each layer. Note that in addition to this expel process ('source' term), !
          !    cumulus updraft can modify layer mean condensate through normal detrainment !
          !    forcing or compensating subsidence.                                         !
          ! 2. Expelled water can be either 'sustaining' or 'precipitating' condensate. By !
          !    definition, 'suataining condensate' will remain in the layer where it was   !
          !    formed, while 'precipitating condensate' will fall across the base of the   !
          !    layer where it was formed.                                                  !
          ! 3. All precipitating condensates are assumed to fall into the release layer or !
          !    ground as soon as it was formed without being evaporated during the falling !
          !    process down to the desinated layer ( either release layer of surface ).    !
          ! ------------------------------------------------------------------------------ !

          ! ------------------------------------------------------------------------- !     
          ! 'dwten(k)','diten(k)' : Production rate of condensate  within the layer k !
          !      [ kg/kg/s ]        by the expels of condensate from cumulus updraft. !
          ! It is important to note that in terms of moisture tendency equation, this !
          ! is a 'source' term of enviromental 'qt'.  More importantly,  these source !
          ! are already counted in the turbulent heat and moisture fluxes we computed !
          ! until now, assuming all the expelled condensate remain in the layer where ! 
          ! it was formed. Thus, in calculation of 'qtten' and 'slten' below, we MUST !
          ! NOT add or subtract these terms explicitly in order not to double or miss !
          ! count, unless some expelled condensates fall down out of the layer.  Note !
          ! this falling-down process ( i.e., precipitation process ) and  associated !
          ! 'qtten' and 'slten' and production of surface precipitation flux  will be !
          ! treated later in 'zm_conv_evap' in 'convect_shallow_tend' subroutine.     ! 
          ! In below, we are converting expelled cloud condensate into correct unit.  !
          ! I found that below use of '0.5 * (umf(k-1) + umf(k))' causes conservation !
          ! errors at some columns in global simulation. So, I returned to originals. !
          ! This will cause no precipitation flux at 'kpen' layer since umf(kpen)=0.  !
          ! ------------------------------------------------------------------------- !

          dwten(k) = dwten(k) * 0.5_r8 * ( umf(k-1) + umf(k) ) * g / dp0(k) ! [ kg/kg/s ]
          diten(k) = diten(k) * 0.5_r8 * ( umf(k-1) + umf(k) ) * g / dp0(k) ! [ kg/kg/s ]  

          ! dwten(k) = dwten(k) * umf(k) * g / dp0(k) ! [ kg/kg/s ]
          ! diten(k) = diten(k) * umf(k) * g / dp0(k) ! [ kg/kg/s ]

          ! --------------------------------------------------------------------------- !
          ! 'qrten(k)','qsten(k)' : Production rate of rain and snow within the layer k !
          !     [ kg/kg/s ]         by cumulus expels of condensates to the environment.!         
          ! This will be falled-out of the layer where it was formed and will be dumped !
          ! dumped into the release layer assuming that there is no evaporative cooling !
          ! while precipitable condensate moves to the relaes level. This is reasonable ! 
          ! assumtion if cumulus is purely vertical and so the path along which precita !
          ! ble condensate falls is fully saturared. This 're-allocation' process of    !
          ! precipitable condensate into the release layer is fully described in this   !
          ! convection scheme. After that, the dumped water into the release layer will !
          ! falling down across the base of release layer ( or LCL, if  exact treatment ! 
          ! is required ) and will be allowed to be evaporated in layers below  release !
          ! layer, and finally non-zero surface precipitation flux will be calculated.  !
          ! This latter process will be separately treated 'zm_conv_evap' routine.      !
          ! --------------------------------------------------------------------------- !

          qrten(k) = frc_rasn * dwten(k)
          qsten(k) = frc_rasn * diten(k) 
 
          ! ----------------------------------------------------------------------- !         
          ! 'rainflx','snowflx' : Cumulative rain and snow flux integrated from the ! 
          !     [ kg/m2/s ]       release leyer to the 'kpen' layer. Note that even !
          ! though wtw(kpen) < 0 (and umf(kpen) = 0) at the top interface of 'kpen' !
          ! layer, 'dwten(kpen)' and diten(kpen)  were calculated after calculating !
          ! explicit cloud top height. Thus below calculation of precipitation flux !
          ! is correct. Note that  precipitating condensates are formed only in the !
          ! layers from 'krel' to 'kpen', including the two layers.                 !
          ! ----------------------------------------------------------------------- !

          rainflx  = rainflx + qrten(k) * dp0(k) / g
          snowflx  = snowflx + qsten(k) * dp0(k) / g

          ! -------------------------------------------------------------------- !
          ! Correct understanding of the meaning of 'qc' is critically important !
          ! because this is directly used in the other parts of subroutines, for !
          ! checking of moisture conservations, and stratiform microhpysics, etc.!
          ! 'qc' is production rate of suspendable condensate within the layer k !
          !      by the sink ( or expel ) of cumulus 'suspended condensate' into ! 
          ! the environment. In terms of moisture tendency equation, this is one !
          ! source of environmental qt ( or sink of cumulus qt ). In other parts !
          ! of CAM, 'qc' is passed by the variable name 'dlf' in stratiform_tend !
          ! ( 'qc' in 'stratiform_tend' actually includes contribution both from !
          ! 'deep' & 'shallow' convection schemes). In other subroutine, this is !
          ! called 'reserved liquid water tendency.                              ! 
          ! Note that 'qrten(k)+qc_l = dwten(k)', 'qsten(k)+qc_i = diten(k)'.    !
          ! -------------------------------------------------------------------- !
 
          qc_l(k)  = ( 1. - frc_rasn ) * dwten(k) ! [ kg/kg/s ]
          qc_i(k)  = ( 1._r8 - frc_rasn ) * diten(k) ! [ kg/kg/s ]
          qc(k) = qc_l(k) + qc_i(k)   

          ! ------------------------------------------------------------------------ !
          ! 'slten(k)','qtten(k)'                                                    !
          !  Note that 'slflx(k)' and 'qtflx(k)' we have calculated already included !
          !  all the contributions of (1) expels of condensate (dwten(k), diten(k)), !
          !  (2) mass detrainment ( delta * umf * ( qtu - qt ) ), & (3) compensating !
          !  subsidence ( M * dqt / dz ). Thus 'slflx(k)' and 'qtflx(k)' we computed ! 
          !  is a hybrid turbulent flux containing one part of 'source' term - expel !
          !  of condensate. In order to calculate 'slten' and 'qtten', we should add !
          !  additional 'source' term, if any. If the expelled condensate falls down !
          !  across the base of the layer, it will be another sink (negative source) !
          !  term.  Note also that we included frictional heating terms in the below !
          !  calculation of 'slten'.                                                 !
          ! ------------------------------------------------------------------------ !
                   
          slten(k) = ( slflx(km1) - slflx(k) ) * g / dp0(k)
          if( k.eq.1 ) then
              slten(k) = slten(k) - g / 4 / dp0(k) * (                                &
                                    uflx(k)*(uf(k+1) - uf(k) + u0(k+1) - u0(k)) +     & 
                                    vflx(k)*(vf(k+1) - vf(k) + v0(k+1) - v0(k)))
          elseif( k.ge.2 .and. k.le.kpen-1 ) then
              slten(k) = slten(k) - g / 4 / dp0(k) * (                                &
                                    uflx(k)*(uf(k+1) - uf(k) + u0(k+1) - u0(k)) +     &
                                    uflx(k-1)*(uf(k) - uf(k-1) + u0(k) - u0(k-1)) +   &
                                    vflx(k)*(vf(k+1) - vf(k) + v0(k+1) - v0(k)) +     &
                                    vflx(k-1)*(vf(k) - vf(k-1) + v0(k) - v0(k-1)))
          elseif( k.eq.kpen ) then
              slten(k) = slten(k) - g / 4 / dp0(k) * (                                &
                                    uflx(k-1)*(uf(k) - uf(k-1) + u0(k) - u0(k-1)) +   &
                                    vflx(k-1)*(vf(k) - vf(k-1) + v0(k) - v0(k-1)))
          endif
          qtten(k) = ( qtflx(km1) - qtflx(k) ) * g / dp0(k)

          ! ------------------------------------------------------------------------- !
          ! 'qlten(k)','qiten(k)','qvten(k)'                                          !
          !  If I want to include 'qlten_sink(k)' and 'qiten_sink(k)' computed above  !
          !  I can simply add this term in the below equation of qlten(k) & qiten(k). !
          !  If either of 'ql' or 'qi' at the next time step is negative, force them  !
          !  to be positive, and calculate qvten(k) accordingly,  in such a way that  !
          !  it conserves total moisture. If qv(k) at the next time step is negative, !
          !  set a minimum value on qv0(k) at the next time step and modify qlten(k)  !
          !  and qiten(k), such that total moisture is conserved.                     !
          ! ------------------------------------------------------------------------- !
           
          qlten(k) = dwten(k) + ( qtten(k) - dwten(k) - diten(k) ) * ( ql0(k) / qt0(k) )
          qiten(k) = diten(k) + ( qtten(k) - dwten(k) - diten(k) ) * ( qi0(k) / qt0(k) )
          qvten(k) = qtten(k) - qlten(k) - qiten(k)

          ! ------------------------------------ !
          !'sten(k) : Dry static energy tendency !
          ! ------------------------------------ !       
 
          sten(k) = slten(k) + xlv * qlten(k) + xls * qiten(k)

          ! -------------------------------------------------------------------------- !
          ! 'rliq' : Verticall-integrated 'suspended cloud condensate'                 !
          !  [m/s]   This is so called 'reserved liquid water'  in other subroutines   ! 
          ! of CAM3, since the contribution of this term should not be included into   !
          ! the tendency of each layer or surface flux (precip)  within this cumulus   !
          ! scheme. The adding of this term to the layer tendency will be done inthe   !
          ! 'stratiform_tend', just after performing sediment process there.           !
          ! The main problem of these rather going-back-and-forth and stupid-seeming   ! 
          ! approach is that the sediment process of suspendened condensate will not   !
          ! be treated at all in the 'stratiform_tend'.                                !
          ! Note that 'precip' [m/s] is vertically-integrated total 'rain+snow' formed !
          ! from the cumulus updraft. Important : in the below, 1000 is rhoh2o ( water !
          ! density ) [ kg/m^3 ] used for unit conversion from [ kg/m^2/s ] to [ m/s ] !
          ! for use in stratiform.F90.                                                 !
          ! -------------------------------------------------------------------------- ! 

          rliq    =  rliq    + qc(k) * dp0(k) / g / 1000._r8    ! [ m/s ]

       end do

          precip  =  rainflx + snowflx                       ! [ kg/m2/s ]
          snow    =  snowflx                                 ! [ kg/m2/s ] 

       ! ---------------------------------------------------------------- !
       ! Now treats the 'evaporation' and 'melting' of rain ( qrten ) and ! 
       ! snow ( qsten ) during falling process. Below algorithms are from !
       ! 'zm_conv_evap' but with some modification, which allows separate !
       ! treatment of 'rain' and 'snow' condensates. Note that I included !
       ! the evaporation dynamics into the convection scheme for complete !
       ! development of cumulus scheme especially in association with the ! 
       ! implicit CIN closure. In compatible with this internal treatment !
       ! of evaporation, I should modify 'convect_shallow',  in such that !
       ! 'zm_conv_evap' is not performed when I choose UW PBL-Cu schemes. !                                          
       ! ---------------------------------------------------------------- !

       evpint_rain = 0._r8 
       evpint_snow = 0._r8
       flxrain(0:mkx) = 0._r8
       flxsnow(0:mkx) = 0._r8
       ntraprd(:mkx) = 0._r8
       ntsnprd(:mkx) = 0._r8

       do k = mkx, 1, -1  ! 'k' is a layer index : 'mkx'('1') is the top ('bottom') layer
          
          ! ----------------------------------------------------------------------------- !
          ! flxsntm [kg/m2/s] : Downward snow flux at the top of each layer after melting.! 
          ! snowmlt [kg/kg/s] : Snow melting tendency.                                    !
          ! Below allows melting of snow when it goes down into the warm layer below.     !
          ! ----------------------------------------------------------------------------- !

          if( t0(k).gt.273.16_r8 ) then
              snowmlt = flxsnow(k) * g / dp0(k) 
          else
              snowmlt = 0._r8
          endif

          ! ----------------------------------------------------------------- !
          ! Evaporation rate of 'rain' and 'snow' in the layer k, [ kg/kg/s ] !
          ! where 'rain' and 'snow' are coming down from the upper layers.    !
          ! I used the same evaporative efficiency both for 'rain' and 'snow'.!
          ! Note that evaporation is not allowed in the layers 'k >= krel' by !
          ! assuming that inside of cumulus cloud, across which precipitation !
          ! is falling down, is fully saturated.                              !
          ! The asumptions in association with the 'evplimit_rain(snow)' are  !
          !   1. Do not allow evaporation to supersate the layer              !
          !   2. Do not evaporate more than the flux falling into the layer   !
          !   3. Total evaporation cannot exceed the input total surface flux !
          ! ----------------------------------------------------------------- !

          status = qsat(t0(k),p0(k),es(1),qs(1),gam(1), 1)          
          subsat = max( ( 1._r8 - qv0(k)/qs(1) ), 0._r8 )
          if( noevap_krelkpen ) then
              if( k.ge.krel ) subsat = 0._r8
          endif
        ! Sungsu
        ! Potentially 'evprain' can contain 'snowmlt' as below as in the unified scheme
        ! evprain  = kevp*(1._r8-cldfrct(k))*subsat*sqrt(flxrain(k)+snowmlt*dp0(k)/g) ! Alternative
          evprain  = kevp*(1._r8-cldfrct(k))*subsat*sqrt(flxrain(k))
          evpsnow  = kevp*(1._r8-cldfrct(k))*subsat*sqrt(max(flxsnow(k)-snowmlt*dp0(k)/g,0._r8))

          evplimit_rain = max( 0._r8, ( qs(1) - qv0(k) ) / dt ) 
          evplimit_rain = min( evplimit_rain, flxrain(k) * g /dp0(k) )
          evplimit_rain = min( evplimit_rain, ( rainflx - evpint_rain ) * g / dp0(k) )
          evprain = min( evplimit_rain, evprain )

          evplimit_snow = max( 0._r8, ( qs(1) - qv0(k) ) / dt ) 
          evplimit_snow = min( evplimit_snow, max(flxsnow(k)-snowmlt*dp0(k)/g ,0._r8) * g /dp0(k) )
          evplimit_snow = min( evplimit_snow, ( snowflx - evpint_snow ) * g / dp0(k) )
          evpsnow = min( evplimit_snow, evpsnow )

          ! ------------------------------------------------------------- !
          ! Vertically-integrated evaporative fluxes of 'rain' and 'snow' !
          ! ------------------------------------------------------------- !

          evpint_rain = evpint_rain + evprain * dp0(k) / g
          evpint_snow = evpint_snow + evpsnow * dp0(k) / g

          ! -------------------------------------------------------------- !
          ! Net 'rain' and 'snow' production rate in the layer [ kg/kg/s ] !
          ! -------------------------------------------------------------- !         

          ntraprd(k) = qrten(k) - evprain + snowmlt
          ntsnprd(k) = qsten(k) - evpsnow - snowmlt
         
          ! -------------------------------------------------------------------------------- !
          ! Downward fluxes of 'rain' and 'snow' fluxes at the base of the layer [ kg/m2/s ] !
          ! Note that layer index increases with height.                                     !
          ! -------------------------------------------------------------------------------- !

          flxrain(k-1) = flxrain(k) + ntraprd(k) * dp0(k) / g
          flxsnow(k-1) = flxsnow(k) + ntsnprd(k) * dp0(k) / g
          flxrain(k-1) = max( flxrain(k-1), 0._r8 )
          if(flxrain(k-1).eq.0._r8) ntraprd(k) = -flxrain(k) * g / dp0(k)
          flxsnow(k-1) = max( flxsnow(k-1), 0._r8 )         
          if(flxsnow(k-1).eq.0._r8) ntsnprd(k) = -flxsnow(k) * g / dp0(k)

          ! ---------------------------------- !
          ! Calculate thermodynamic tendencies !
          ! --------------------------------------------------------------------------- !
          ! Note that equivalently, we can write tendency formula of 'sten' and 'slten' !
          ! by 'sten(k)  = sten(k) - xlv*evprain  - xls*evpsnow - (xls-xlv)*snowmlt' &  !
          !    'slten(k) = sten(k) - xlv*qlten(k) - xls*qiten(k)'.                      !
          ! The above formula is equivalent to the below formula. However below formula !
          ! is preferred since we have already imposed explicit constraint on 'ntraprd' !
          ! and 'ntsnprd' in case that flxrain(k-1) < 0 & flxsnow(k-1) < 0._r8             !
          ! Note : In future, I can elborate the limiting of 'qlten','qvten','qiten'    !
          !        such that that energy and moisture conservation error is completely  !
          !        suppressed.                                                          !
          ! --------------------------------------------------------------------------- !

          qlten(k) = max( qlten(k) - qrten(k), ( 0.e-16_r8 - ql0(k) ) / dt )
          qiten(k) = max( qiten(k) - qsten(k), ( 0.e-16_r8 - qi0(k) ) / dt )
          qvten(k) = max( qvten(k) + evprain + evpsnow, ( 2.e-12_r8 - qv0(k) ) / dt )
          qtten(k) = qlten(k) + qiten(k) + qvten(k)
          if(qvten(k).eq.((2.e-12_r8-qv0(k))/dt).or.qlten(k).eq.(-ql0(k)/dt) &
                                             .or.qiten(k).eq.(-qi0(k)/dt)) then
             limit_negcon(i) = 1._r8
          end if
          slten(k) = slten(k) + xlv * ntraprd(k) + xls * ntsnprd(k)         
          sten(k)  = slten(k) + xlv * qlten(k)   + xls * qiten(k)

       end do

       ! ------------------------------------------------------------- !
       ! Calculate final surface flux of precipitation, rain, and snow !
       ! Convert unit to [m/s] for use in 'check_energy_chng'.         !  
       ! ------------------------------------------------------------- !

       precip  = ( flxrain(0) + flxsnow(0) ) / 1000._r8
       snow    =   flxsnow(0) / 1000._r8       

       ! --------------------------------------------------------------------------- !
       ! Until now, all the calculations are done completely in this shallow cumulus !
       ! scheme. If you want to use this cumulus scheme other than CAM3, then do not !
       ! perform below block. However, for compatible use with the other subroutines !
       ! in CAM3, I should subtract the effect of 'qc(k)' ('rliq') from the tendency !
       ! equation in each layer, since this effect will be separately added later in !
       ! in 'stratiform_tend' just after performing sediment process there. In order !
       ! to be consistent with 'stratiform_tend', just subtract qc(k)  from tendency !
       ! equation of each layer, but do not add it to the 'precip'. Apprently,  this !
       ! will violate energy and moisture conservations.    However, when performing !
       ! conservation check in 'tphysbc.F90' just after 'convect_shallow_tend',   we !
       ! will add 'qc(k)' ( rliq ) to the surface flux term just for the purpose  of !
       ! passing the energy-moisture conservation check. Explicit adding-back of 'qc'!
       ! to the individual layer tendency equation will be done in 'stratiform_tend' !
       ! after performing sediment process there. Simply speaking, in 'tphysbc' just !
       ! after 'convect_shallow_tend', we will dump 'rliq' into surface as a  'rain' !
       ! in order to satisfy energy and moisture conservation, and  in the following !
       ! 'stratiform_tend', we will restore it back to 'qlten(k)' ( 'ice' will go to !  
       ! 'water' there) from surface precipitation. This is a funny but conceptually !
       ! entertaining procedure. One concern I have for this complex process is that !
       ! output-writed stratiform precipitation amount will be underestimated due to !
       ! arbitrary subtracting of 'rliq' in stratiform_tend, where                   !
       ! ' prec_str = prec_sed + prec_pcw - rliq' and 'rliq' is not real but fake.   ! 
       ! However, as shown in 'srfxfer.F90', large scale precipitation amount (PRECL)!
       ! that is writed-output is corrected written since in 'srfxfer.F90',  PRECL = !
       ! 'prec_sed + prec_pcw', without including 'rliq'. So current code is correct.!
       ! Note also in 'srfxfer.F90', convective precipitation amount is 'PRECC =     ! 
       ! prec_zmc(i) + prec_cmf(i)' which is also correct.                           !
       ! --------------------------------------------------------------------------- !

       do k = 1, kpen       
          qtten(k) = qtten(k) - qc(k)
          qlten(k) = qlten(k) - qc_l(k)
          qiten(k) = qiten(k) - qc_i(k)
          slten(k) = slten(k) + ( xlv * qc_l(k) + xls * qc_i(k) )
          ! ---------------------------------------------------------------------- !
          ! Since all reserved condensates will be treated  as liquid water in the !
          ! 'check_energy_chng' & 'stratiform_tend' without an explicit conversion !
          ! algorithm, I should consider explicitly the energy conversions between !
          ! 'ice' and 'liquid' - i.e., I should convert 'ice' to 'liquid'  and the !
          ! necessary energy for this conversion should be subtracted from 'sten'. ! 
          ! Without this conversion here, energy conservation error come out. Note !
          ! that there should be no change of 'qvten(k)'.                          !
          ! ---------------------------------------------------------------------- !
          sten(k)  = sten(k)  - ( xls - xlv ) * qc_i(k)
       end do
       
       ! ---------------------------------------------------------------- !
       ! Cumpute default diagnostic outputs                               !
       ! Note that since 'qtu(krel-1:kpen-1)' & 'thlu(krel-1:kpen-1)' has !
       ! been adjusted after detraining cloud condensate into environment ! 
       ! during cumulus updraft motion,  below calculations will  exactly !
       ! reproduce in-cloud properties as shown in the output analysis.   !
       ! ---------------------------------------------------------------- ! 
 
       call conden(prel,thlu(krel-1),qtu(krel-1),thj,qvj,qlj,qij,qse,id_check,qsat)
       if(id_check.eq.1) then
          exit_conden(i) = 1._r8
          id_exit = .true.
          go to 333
       end if
       qcubelow = qlj + qij
       qlubelow = qlj       
       qiubelow = qij       
       rcwp = 0._r8
       rlwp = 0._r8
       riwp = 0._r8

       ! --------------------------------------------------------------------- !
       ! In the below calculations, I explicitly considered cloud base ( LCL ) !
       ! and cloud top height ( ps0(kpen-1) + ppen )                           !
       ! ----------------------------------------------------------------------! 
       do k = krel, kpen ! This is a layer index
          ! ------------------------------------------------------------------ ! 
          ! Calculate cumulus condensate at the upper interface of each layer. !
          ! Note 'ppen < 0' and at 'k=kpen' layer, I used 'thlu_top'&'qtu_top' !
          ! which explicitly considered zero or non-zero 'fer(kpen)'.          !
          ! ------------------------------------------------------------------ ! 
          if( k.eq.kpen ) then 
              call conden(ps0(k-1)+ppen,thlu_top,qtu_top,thj,qvj,qlj,qij,qse,id_check,qsat)
          else
              call conden(ps0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check,qsat)
          endif
          if(id_check.eq.1) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
          end if
          ! ---------------------------------------------------------------- !
          ! Calculate in-cloud mean LWC ( qlu(k) ), IWC ( qiu(k) ),  & layer !
          ! mean cumulus fraction ( cufrc(k) ),  vertically-integrated layer !
          ! mean LWP and IWP. Expel some of in-cloud condensate at the upper !
          ! interface if it is largr than criqc. Note cumulus cloud fraction !
          ! is assumed to be twice of core updraft fractional area. Thus LWP !
          ! and IWP will be twice of actual value coming from our scheme.    !
          ! ---------------------------------------------------------------- !
          qcu(k) = 0.5_r8 * ( qcubelow + qlj + qij )
          qlu(k) = 0.5_r8 * ( qlubelow + qlj )
          qiu(k) = 0.5_r8 * ( qiubelow + qij )
          cufrc(k) = ( ufrc(k-1) + ufrc(k) )
          if( k.eq.krel ) then
              cufrc(k) = ( ufrclcl + ufrc(k) )*( prel - ps0(k) )/(ps0(k-1)-ps0(k))
          else if( k.eq.kpen ) then
              cufrc(k) = ( ufrc(k-1) +    0._r8 )*( -ppen )        /(ps0(k-1)-ps0(k))
              if( (qlj + qij) .gt. criqc ) then           
                   qcu(k) = 0.5_r8 * ( qcubelow + criqc )
                   qlu(k) = 0.5_r8 * ( qlubelow + criqc * qlj / ( qlj + qij ) )
                   qiu(k) = 0.5_r8 * ( qiubelow + criqc * qij / ( qlj + qij ) )
              endif
          endif  
          rcwp = rcwp + ( qlu(k) + qiu(k) ) * ( ps0(k-1) - ps0(k) ) / g * cufrc(k)
          rlwp = rlwp +   qlu(k)            * ( ps0(k-1) - ps0(k) ) / g * cufrc(k)
          riwp = riwp +   qiu(k)            * ( ps0(k-1) - ps0(k) ) / g * cufrc(k)
          qcubelow = qlj + qij
          qlubelow = qlj
          qiubelow = qij
       end do
       ! ------------------------------------ !      
       ! Cloud top and base interface indices !
       ! ------------------------------------ !
       cnt = real(kpen, r8)
       cnb = real(krel - 1, r8)

       ! ------------------------------------------------------------------------- !
       ! End of formal calculation. Below blocks are for implicit CIN calculations ! 
       ! with re-initialization and save variables at iter_cin = 1._r8                !
       ! ------------------------------------------------------------------------- !
       
       ! --------------------------------------------------------------- !
       ! Adjust the original input profiles for implicit CIN calculation !
       ! --------------------------------------------------------------- !

       if(iter .ne. iter_cin) then 

          ! ------------------------------------------------------------------- !
          ! Save the output from "iter_cin = 1"                                 !
          ! These output will be writed-out if "iter_cin = 1" was not performed !
          ! for some reasons.                                                   !
          ! ------------------------------------------------------------------- !

          qv0_s(:mkx)           = qv0(:mkx) + qvten(:mkx) * dt
          ql0_s(:mkx)           = ql0(:mkx) + qlten(:mkx) * dt
          qi0_s(:mkx)           = qi0(:mkx) + qiten(:mkx) * dt
          s0_s(:mkx)            = s0(:mkx) + sten(:mkx) * dt 
          u0_s(:mkx)            = u0(:mkx) + uten(:mkx) * dt
          v0_s(:mkx)            = v0(:mkx) + vten(:mkx) * dt 
          qt0_s(:mkx)           = qv0_s(:mkx) + ql0_s(:mkx) + qi0_s(:mkx)
          t0_s(:mkx)            = t0(:mkx) + sten(:mkx) * dt / cp

          umf_s(0:mkx)          = umf(0:mkx)
          qvten_s(:mkx)         = qvten(:mkx)
          qlten_s(:mkx)         = qlten(:mkx)  
          qiten_s(:mkx)         = qiten(:mkx)
          sten_s(:mkx)          = sten(:mkx)
          uten_s(:mkx)          = uten(:mkx)  
          vten_s(:mkx)          = vten(:mkx)
          qrten_s(:mkx)         = qrten(:mkx)
          qsten_s(:mkx)         = qsten(:mkx)  
          precip_s              = precip
          snow_s                = snow
          cush_s                = cush
          cufrc_s(:mkx)         = cufrc(:mkx)  
          slflx_s(0:mkx)        = slflx(0:mkx)  
          qtflx_s(0:mkx)        = qtflx(0:mkx)  
          qcu_s(:mkx)           = qcu(:mkx)  
          qlu_s(:mkx)           = qlu(:mkx)  
          qiu_s(:mkx)           = qiu(:mkx)  
          fer_s(:mkx)           = fer(:mkx)  
          fdr_s(:mkx)           = fdr(:mkx)  
          cin_s                 = cin
          cinlcl_s              = cinlcl
          cbmf_s                = cbmf
          rliq_s                = rliq
          qc_s(:mkx)            = qc(:mkx)
          cnt_s                 = cnt
          cnb_s                 = cnb
          qtten_s(:mkx)         = qtten(:mkx)
          slten_s(:mkx)         = slten(:mkx)
          ufrc_s(0:mkx)         = ufrc(0:mkx) 

          uflx_s(0:mkx)         = uflx(0:mkx)  
          vflx_s(0:mkx)         = vflx(0:mkx)  
          
          ufrcinvbase_s         = ufrcinvbase
          ufrclcl_s             = ufrclcl 
          winvbase_s            = winvbase
          wlcl_s                = wlcl
          plcl_s                = plcl
          pinv_s                = ps0(kinv-1)
          plfc_s                = plfc        
          pbup_s                = ps0(kbup)
          ppen_s                = ps0(kpen-1)+ppen        
          qtsrc_s               = qtsrc
          thlsrc_s              = thlsrc
          thvlsrc_s             = thvlsrc
          emfkbup_s             = emf(kbup)
          cbmflimit_s           = cbmflimit
          tkeavg_s              = tkeavg
          zinv_s                = zs0(kinv-1)
          rcwp_s                = rcwp
          rlwp_s                = rlwp
          riwp_s                = riwp

          wu_s(0:mkx)           = wu(0:mkx)
          qtu_s(0:mkx)          = qtu(0:mkx)
          thlu_s(0:mkx)         = thlu(0:mkx)
          thvu_s(0:mkx)         = thvu(0:mkx)
          uu_s(0:mkx)           = uu(0:mkx)
          vu_s(0:mkx)           = vu(0:mkx)
          qtu_emf_s(0:mkx)      = qtu_emf(0:mkx)
          thlu_emf_s(0:mkx)     = thlu_emf(0:mkx)
          uu_emf_s(0:mkx)       = uu_emf(0:mkx)
          vu_emf_s(0:mkx)       = vu_emf(0:mkx)
          uemf_s(0:mkx)         = uemf(0:mkx)

          dwten_s(:mkx)         = dwten(:mkx)
          diten_s(:mkx)         = diten(:mkx)
          flxrain_s(0:mkx)      = flxrain(0:mkx)
          flxsnow_s(0:mkx)      = flxsnow(0:mkx)
          ntraprd_s(:mkx)       = ntraprd(:mkx)
          ntsnprd_s(:mkx)       = ntsnprd(:mkx)

          excessu_arr_s(:mkx)   = excessu_arr(:mkx)
          excess0_arr_s(:mkx)   = excess0_arr(:mkx)
          xc_arr_s(:mkx)        = xc_arr(:mkx)
          aquad_arr_s(:mkx)     = aquad_arr(:mkx)
          bquad_arr_s(:mkx)     = bquad_arr(:mkx)
          cquad_arr_s(:mkx)     = cquad_arr(:mkx)
          bogbot_arr_s(:mkx)    = bogbot_arr(:mkx)
          bogtop_arr_s(:mkx)    = bogtop_arr(:mkx)

          ! ----------------------------------------------------------------------------- ! 
          ! Recalculate environmental variables for new cin calculation at "iter_cin = 2" ! 
          ! using the updated state variables. Perform only for variables necessary  for  !
          ! the new cin calculation.                                                      !
          ! ----------------------------------------------------------------------------- !
          
          qv0(:mkx) = qv0_s(:mkx)
          ql0(:mkx) = ql0_s(:mkx)
          qi0(:mkx) = qi0_s(:mkx)
          s0(:mkx)  = s0_s(:mkx)
          t0(:mkx)  = t0_s(:mkx)
      
          qt0(:mkx)   = (qv0(:mkx) + ql0(:mkx) + qi0(:mkx))
          thl0(:mkx)  = (t0(:mkx) - xlv*ql0(:mkx)/cp - xls*qi0(:mkx)/cp)/exn0(:mkx)
          thvl0(:mkx) = (1._r8 + zvir*qt0(:mkx))*thl0(:mkx)

          ssthl0 = slope(mkx,thl0,p0) ! Dimension of ssthl0(:mkx) is implicit
          ssqt0  = slope(mkx,qt0 ,p0)
          ssu0   = slope(mkx,u0  ,p0)
          ssv0   = slope(mkx,v0  ,p0)

          do k = 1, mkx

          thl0bot = thl0(k) + ssthl0(k)*(ps0(k-1) - p0(k))
          qt0bot  = qt0(k)  + ssqt0(k) *(ps0(k-1) - p0(k))
          call conden(ps0(k-1),thl0bot,qt0bot,thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
          end if
          thv0bot(k)  = thj*(1._r8 + zvir*qvj - qlj - qij)
          thvl0bot(k) = thl0bot*(1._r8 + zvir*qt0bot)
          
          thl0top     = thl0(k) + ssthl0(k)*(ps0(k) - p0(k))
          qt0top      =  qt0(k) + ssqt0(k) *(ps0(k) - p0(k))
          call conden(ps0(k),thl0top,qt0top,thj,qvj,qlj,qij,qse,id_check,qsat)
          if(id_check.eq.1) then
             exit_conden(i) = 1._r8
             id_exit = .true.
             go to 333
          end if
          thv0top(k)  = thj*(1._r8 + zvir*qvj - qlj - qij)
          thvl0top(k) = thl0top*(1._r8 + zvir*qt0top)

          end do

       endif               ! End of 'if(iter .ne. iter_cin)' if sentence. 

     end do                ! End of implicit CIN loop (cin_iter)      

     ! ----------------------- !
     ! Update Output Variables !
     ! ----------------------- !

     umf_out(i,0:mkx)         = umf(0:mkx)
     slflx_out(i,0:mkx)       = slflx(0:mkx)
     qtflx_out(i,0:mkx)       = qtflx(0:mkx)
     qvten_out(i,:mkx)        = qvten(:mkx)
     qlten_out(i,:mkx)        = qlten(:mkx)
     qiten_out(i,:mkx)        = qiten(:mkx)
     sten_out(i,:mkx)         = sten(:mkx)
     uten_out(i,:mkx)         = uten(:mkx)
     vten_out(i,:mkx)         = vten(:mkx)
     qrten_out(i,:mkx)        = qrten(:mkx)
     qsten_out(i,:mkx)        = qsten(:mkx)
     precip_out(i)            = precip
     snow_out(i)              = snow
     cufrc_out(i,:mkx)        = cufrc(:mkx)
     qcu_out(i,:mkx)          = qcu(:mkx)
     qlu_out(i,:mkx)          = qlu(:mkx)
     qiu_out(i,:mkx)          = qiu(:mkx)
     cush_inout(i)            = cush
     cbmf_out(i)              = cbmf
     rliq_out(i)              = rliq
     qc_out(i,:mkx)           = qc(:mkx)
     cnt_out(i)               = cnt
     cnb_out(i)               = cnb
  
     ! ------------------------------------------------- !
     ! Below are specific diagnostic output for detailed !
     ! analysis of cumulus scheme                        !
     ! ------------------------------------------------- !

     fer_out(i,mkx:1:-1)      = fer(:mkx)  
     fdr_out(i,mkx:1:-1)      = fdr(:mkx)  
     cinh_out(i)              = cin
     cinlclh_out(i)           = cinlcl
     qtten_out(i,mkx:1:-1)    = qtten(:mkx)
     slten_out(i,mkx:1:-1)    = slten(:mkx)
     ufrc_out(i,mkx:0:-1)     = ufrc(0:mkx)
     uflx_out(i,mkx:0:-1)     = uflx(0:mkx)  
     vflx_out(i,mkx:0:-1)     = vflx(0:mkx)  
     
     ufrcinvbase_out(i)       = ufrcinvbase
     ufrclcl_out(i)           = ufrclcl 
     winvbase_out(i)          = winvbase
     wlcl_out(i)              = wlcl
     plcl_out(i)              = plcl
     pinv_out(i)              = ps0(kinv-1)
     plfc_out(i)              = plfc    
     pbup_out(i)              = ps0(kbup)        
     ppen_out(i)              = ps0(kpen-1)+ppen            
     qtsrc_out(i)             = qtsrc
     thlsrc_out(i)            = thlsrc
     thvlsrc_out(i)           = thvlsrc
     emfkbup_out(i)           = emf(kbup)
     cbmflimit_out(i)         = cbmflimit
     tkeavg_out(i)            = tkeavg
     zinv_out(i)              = zs0(kinv-1)
     rcwp_out(i)              = rcwp
     rlwp_out(i)              = rlwp
     riwp_out(i)              = riwp

     wu_out(i,mkx:0:-1)       = wu(0:mkx)
     qtu_out(i,mkx:0:-1)      = qtu(0:mkx)
     thlu_out(i,mkx:0:-1)     = thlu(0:mkx)
     thvu_out(i,mkx:0:-1)     = thvu(0:mkx)
     uu_out(i,mkx:0:-1)       = uu(0:mkx)
     vu_out(i,mkx:0:-1)       = vu(0:mkx)
     qtu_emf_out(i,mkx:0:-1)  = qtu_emf(0:mkx)
     thlu_emf_out(i,mkx:0:-1) = thlu_emf(0:mkx)
     uu_emf_out(i,mkx:0:-1)   = uu_emf(0:mkx)
     vu_emf_out(i,mkx:0:-1)   = vu_emf(0:mkx)
     uemf_out(i,mkx:0:-1)     = uemf(0:mkx)

     dwten_out(i,mkx:1:-1)    = dwten(:mkx)
     diten_out(i,mkx:1:-1)    = diten(:mkx)
     flxrain_out(i,mkx:0:-1)  = flxrain(0:mkx)
     flxsnow_out(i,mkx:0:-1)  = flxsnow(0:mkx)
     ntraprd_out(i,mkx:1:-1)  = ntraprd(:mkx)
     ntsnprd_out(i,mkx:1:-1)  = ntsnprd(:mkx)

     excessu_arr_out(i,mkx:1:-1)  = excessu_arr(:mkx)
     excess0_arr_out(i,mkx:1:-1)  = excess0_arr(:mkx)
     xc_arr_out(i,mkx:1:-1)       = xc_arr(:mkx)
     aquad_arr_out(i,mkx:1:-1)    = aquad_arr(:mkx)
     bquad_arr_out(i,mkx:1:-1)    = bquad_arr(:mkx)
     cquad_arr_out(i,mkx:1:-1)    = cquad_arr(:mkx)
     bogbot_arr_out(i,mkx:1:-1)   = bogbot_arr(:mkx)
     bogtop_arr_out(i,mkx:1:-1)   = bogtop_arr(:mkx)

 333 if(id_exit) then ! Exit without cumulus convection

     exit_UWCu(i) = 1._r8

     ! --------------------------------------------------------------------- !
     ! Initialize output variables when cumulus convection was not performed.!
     ! --------------------------------------------------------------------- !
     
     ! ---------------------------------------------------------------- !
     ! Below block, with !! lines in the below are diagnostic  purpose  !
     ! to look at values even when 333 exit mainly for SCAM simulation. !
     ! When running global simulations, it is good to initialize by  0, !
     ! if cumulus convection was not performed.                         !
     ! ---------------------------------------------------------------- !

     !! umf_out(i,0:mkx)              = umf(0:mkx)
     !! cinh_out(i)                   = cin
     !! cinlclh_out(i)                = cinlcl
     !! cbmf_out(i)                   = cbmf
     !! ufrcinvbase_out(i)            = ufrcinvbase
     !! ufrclcl_out(i)                = ufrclcl 
     !! winvbase_out(i)               = winvbase
     !! wlcl_out(i)                   = wlcl
     !! plcl_out(i)                   = plcl
     !! pinv_out(i)                   = ps0(kinv-1)
     !! plfc_out(i)                   = plfc    
     !! qtsrc_out(i)                  = qtsrc
     !! thlsrc_out(i)                 = thlsrc
     !! thvlsrc_out(i)                = thvlsrc
     !! cbmflimit_out(i)              = cbmflimit
     !! tkeavg_out(i)                 = tkeavg
     !! zinv_out(i)                   = zs0(kinv-1)
     !! rcwp_out(i)                   = rcwp
     !! rlwp_out(i)                   = rlwp
     !! riwp_out(i)                   = riwp

     !! wu_out(i,mkx:0:-1)            = wu(0:mkx)
     !! qtu_out(i,mkx:0:-1)           = qtu(0:mkx)
     !! thlu_out(i,mkx:0:-1)          = thlu(0:mkx)
     !! thvu_out(i,mkx:0:-1)          = thvu(0:mkx)
     !! uu_out(i,mkx:0:-1)            = uu(0:mkx)
     !! vu_out(i,mkx:0:-1)            = vu(0:mkx)
     !! qtu_emf_out(i,mkx:0:-1)       = qtu_emf(0:mkx)
     !! thlu_emf_out(i,mkx:0:-1)      = thlu_emf(0:mkx)
     !! uu_emf_out(i,mkx:0:-1)        = uu_emf(0:mkx)
     !! vu_emf_out(i,mkx:0:-1)        = vu_emf(0:mkx)
     !! uemf_out(i,mkx:0:-1)          = uemf(0:mkx)

     !! dwten_out(i,mkx:1:-1)         = dwten(:mkx)
     !! diten_out(i,mkx:1:-1)         = diten(:mkx)
     !! flxrain_out(i,mkx:0:-1)       = flxrain(0:mkx)
     !! flxsnow_out(i,mkx:0:-1)       = flxsnow(0:mkx)
     !! ntraprd_out(i,mkx:1:-1)       = ntraprd(:mkx)
     !! ntsnprd_out(i,mkx:1:-1)       = ntsnprd(:mkx)

     !! excessu_arr_out(i,mkx:1:-1)   = excessu_arr(:mkx)
     !! excess0_arr_out(i,mkx:1:-1)   = excess0_arr(:mkx)
     !! xc_arr_out(i,mkx:1:-1)        = xc_arr(:mkx)
     !! aquad_arr_out(i,mkx:1:-1)     = aquad_arr(:mkx)
     !! bquad_arr_out(i,mkx:1:-1)     = bquad_arr(:mkx)
     !! cquad_arr_out(i,mkx:1:-1)     = cquad_arr(:mkx)
     !! bogbot_arr_out(i,mkx:1:-1)    = bogbot_arr(:mkx)
     !! bogtop_arr_out(i,mkx:1:-1)    = bogtop_arr(:mkx)

     ! ---------------------------------------------------------- !
     ! Below block is normal initialization for global simulation !
     ! ---------------------------------------------------------- !
     
     umf_out(i,0:mkx)              = 0._r8    !!
     slflx_out(i,0:mkx)            = 0._r8
     qtflx_out(i,0:mkx)            = 0._r8
     qvten_out(i,:mkx)             = 0._r8
     qlten_out(i,:mkx)             = 0._r8
     qiten_out(i,:mkx)             = 0._r8
     sten_out(i,:mkx)              = 0._r8
     uten_out(i,:mkx)              = 0._r8
     vten_out(i,:mkx)              = 0._r8
     qrten_out(i,:mkx)             = 0._r8
     qsten_out(i,:mkx)             = 0._r8
     precip_out(i)                 = 0._r8
     snow_out(i)                   = 0._r8
     cufrc_out(i,:mkx)             = 0._r8
     qcu_out(i,:mkx)               = 0._r8
     qlu_out(i,:mkx)               = 0._r8
     qiu_out(i,:mkx)               = 0._r8
     cush_inout(i)                 = -1._r8
     cbmf_out(i)                   = 0._r8    !!
     rliq_out(i)                   = 0._r8
     qc_out(i,:mkx)                = 0._r8
     ! Set CLDTOP to bottom layer index and CLDBOT to top layer index in
     ! cloud free columns.
     cnt_out(i)                    = 1._r8
     cnb_out(i)                    = real(mkx, r8)

     fer_out(i,mkx:1:-1)           = 0._r8  
     fdr_out(i,mkx:1:-1)           = 0._r8  
     cinh_out(i)                   = -1._r8   !!
     cinlclh_out(i)                = -1._r8   !!
     qtten_out(i,mkx:1:-1)         = 0._r8
     slten_out(i,mkx:1:-1)         = 0._r8
     ufrc_out(i,mkx:0:-1)          = 0._r8
     uflx_out(i,mkx:0:-1)          = 0._r8  
     vflx_out(i,mkx:0:-1)          = 0._r8  

     ufrcinvbase_out(i)            = 0._r8    !!
     ufrclcl_out(i)                = 0._r8    !!
     winvbase_out(i)               = 0._r8    !!
     wlcl_out(i)                   = 0._r8    !!
     plcl_out(i)                   = 0._r8    !!
     pinv_out(i)                   = 0._r8    !! 
     plfc_out(i)                   = 0._r8    !! 
     pbup_out(i)                   = 0._r8    
     ppen_out(i)                   = 0._r8    
     qtsrc_out(i)                  = 0._r8    !!
     thlsrc_out(i)                 = 0._r8    !!
     thvlsrc_out(i)                = 0._r8    !!
     emfkbup_out(i)                = 0._r8
     cbmflimit_out(i)              = 0._r8    !!
     tkeavg_out(i)                 = 0._r8    !!
     zinv_out(i)                   = 0._r8    !!
     rcwp_out(i)                   = 0._r8    !!
     rlwp_out(i)                   = 0._r8    !!
     riwp_out(i)                   = 0._r8    !!

     wu_out(i,mkx:0:-1)            = 0._r8    !!
     qtu_out(i,mkx:0:-1)           = 0._r8    !!    
     thlu_out(i,mkx:0:-1)          = 0._r8    !!     
     thvu_out(i,mkx:0:-1)          = 0._r8    !!     
     uu_out(i,mkx:0:-1)            = 0._r8    !!    
     vu_out(i,mkx:0:-1)            = 0._r8    !!    
     qtu_emf_out(i,mkx:0:-1)       = 0._r8    !!     
     thlu_emf_out(i,mkx:0:-1)      = 0._r8    !!     
     uu_emf_out(i,mkx:0:-1)        = 0._r8    !!      
     vu_emf_out(i,mkx:0:-1)        = 0._r8    !!
     uemf_out(i,mkx:0:-1)          = 0._r8    !!
   
     dwten_out(i,mkx:1:-1)         = 0._r8    !!
     diten_out(i,mkx:1:-1)         = 0._r8    !!
     flxrain_out(i,mkx:0:-1)       = 0._r8    !! 
     flxsnow_out(i,mkx:0:-1)       = 0._r8    !!
     ntraprd_out(i,mkx:1:-1)       = 0._r8    !!
     ntsnprd_out(i,mkx:1:-1)       = 0._r8    !!

     excessu_arr_out(i,mkx:1:-1)   = 0._r8    !!
     excess0_arr_out(i,mkx:1:-1)   = 0._r8    !!
     xc_arr_out(i,mkx:1:-1)        = 0._r8    !!
     aquad_arr_out(i,mkx:1:-1)     = 0._r8    !!
     bquad_arr_out(i,mkx:1:-1)     = 0._r8    !!
     cquad_arr_out(i,mkx:1:-1)     = 0._r8    !!
     bogbot_arr_out(i,mkx:1:-1)    = 0._r8    !!
     bogtop_arr_out(i,mkx:1:-1)    = 0._r8    !!

     end if


     end do                  ! end of big i loop for each column.


     ! ---------------------------------------- !
     ! Writing main diagnostic output variables !
     ! ---------------------------------------- !

     call outfld('qtflx_Cu', qtflx_out(:,mkx:0:-1), mix, lchnk) 
     call outfld('slflx_Cu', slflx_out(:,mkx:0:-1), mix, lchnk) 
     call outfld('uflx_Cu', uflx_out, mix, lchnk) 
     call outfld('vflx_Cu', vflx_out, mix, lchnk) 

     call outfld('qtten_Cu', qtten_out, mix, lchnk) 
     call outfld('slten_Cu', slten_out, mix, lchnk) 
     call outfld('uten_Cu', uten_out(:,mkx:1:-1), mix, lchnk) 
     call outfld('vten_Cu', vten_out(:,mkx:1:-1), mix, lchnk) 
     call outfld('qvten_Cu', qvten_out(:,mkx:1:-1), mix, lchnk) 
     call outfld('qlten_Cu', qlten_out(:,mkx:1:-1), mix, lchnk)
     call outfld('qiten_Cu', qiten_out(:,mkx:1:-1), mix, lchnk)  

     call outfld('cbmf_Cu', cbmf_out, mix, lchnk) 
     call outfld('ufrcinvbase_Cu', ufrcinvbase_out, mix, lchnk) 
     call outfld('ufrclcl_Cu', ufrclcl_out, mix, lchnk) 
     call outfld('winvbase_Cu', winvbase_out, mix, lchnk) 
     call outfld('wlcl_Cu', wlcl_out, mix, lchnk) 
     call outfld('plcl_Cu', plcl_out, mix, lchnk) 
     call outfld('pinv_Cu', pinv_out, mix, lchnk) 
     call outfld('plfc_Cu', plfc_out, mix, lchnk) 
     call outfld('pbup_Cu', pbup_out, mix, lchnk) 
     call outfld('ppen_Cu', ppen_out, mix, lchnk) 
     call outfld('qtsrc_Cu', qtsrc_out, mix, lchnk) 
     call outfld('thlsrc_Cu', thlsrc_out, mix, lchnk) 
     call outfld('thvlsrc_Cu', thvlsrc_out, mix, lchnk) 
     call outfld('emfkbup_Cu', emfkbup_out, mix, lchnk)
     call outfld('cin_Cu', cinh_out, mix, lchnk)  
     call outfld('cinlcl_Cu', cinlclh_out, mix, lchnk) 
     call outfld('cbmflimit_Cu', cbmflimit_out, mix, lchnk) 
     call outfld('tkeavg_Cu', tkeavg_out, mix, lchnk)
     call outfld('zinv_Cu', zinv_out, mix, lchnk)  
     call outfld('rcwp_Cu', rcwp_out, mix, lchnk)
     call outfld('rlwp_Cu', rlwp_out, mix, lchnk)
     call outfld('riwp_Cu', riwp_out, mix, lchnk)
     call outfld('tophgt_Cu', cush_inout, mix, lchnk)   

     call outfld('wu_Cu', wu_out, mix, lchnk)
     call outfld('ufrc_Cu', ufrc_out, mix, lchnk)
     call outfld('qtu_Cu', qtu_out, mix, lchnk)
     call outfld('thlu_Cu', thlu_out, mix, lchnk)
     call outfld('thvu_Cu', thvu_out, mix, lchnk)
     call outfld('uu_Cu', uu_out, mix, lchnk)
     call outfld('vu_Cu', vu_out, mix, lchnk)
     call outfld('qtu_emf_Cu', qtu_emf_out, mix, lchnk)
     call outfld('thlu_emf_Cu', thlu_emf_out, mix, lchnk)
     call outfld('uu_emf_Cu', uu_emf_out, mix, lchnk)
     call outfld('vu_emf_Cu', vu_emf_out, mix, lchnk)
     call outfld('umf_Cu', umf_out(:,mkx:0:-1), mix, lchnk)
     call outfld('uemf_Cu', uemf_out, mix, lchnk)
     call outfld('qcu_Cu', qcu_out(:,mkx:1:-1), mix, lchnk)
     call outfld('qlu_Cu', qlu_out(:,mkx:1:-1), mix, lchnk)
     call outfld('qiu_Cu', qiu_out(:,mkx:1:-1), mix, lchnk)
     call outfld('cufrc_Cu', cufrc_out(:,mkx:1:-1), mix, lchnk)  
     call outfld('fer_Cu', fer_out, mix, lchnk)  
     call outfld('fdr_Cu', fdr_out, mix, lchnk)  

     call outfld('dwten_Cu', dwten_out, mix, lchnk)
     call outfld('diten_Cu', diten_out, mix, lchnk)
     call outfld('qrten_Cu', qrten_out(:,mkx:1:-1), mix, lchnk)
     call outfld('qsten_Cu', qsten_out(:,mkx:1:-1), mix, lchnk)
     call outfld('flxrain_Cu', flxrain_out, mix, lchnk)
     call outfld('flxsnow_Cu', flxsnow_out, mix, lchnk)
     call outfld('ntraprd_Cu', ntraprd_out, mix, lchnk)
     call outfld('ntsnprd_Cu', ntsnprd_out, mix, lchnk)

     call outfld('excessu_Cu', excessu_arr_out, mix, lchnk)
     call outfld('excess0_Cu', excess0_arr_out, mix, lchnk)
     call outfld('xc_Cu', xc_arr_out, mix, lchnk)
     call outfld('aquad_Cu', aquad_arr_out, mix, lchnk)
     call outfld('bquad_Cu', bquad_arr_out, mix, lchnk)
     call outfld('cquad_Cu', cquad_arr_out, mix, lchnk)
     call outfld('bogbot_Cu', bogbot_arr_out, mix, lchnk)
     call outfld('bogtop_Cu', bogtop_arr_out, mix, lchnk)

     call outfld('exit_UWCu_Cu', exit_UWCu, mix, lchnk) 
     call outfld('exit_conden_Cu', exit_conden, mix, lchnk) 
     call outfld('exit_klclmkx_Cu', exit_klclmkx, mix, lchnk) 
     call outfld('exit_klfcmkx_Cu', exit_klfcmkx, mix, lchnk) 
     call outfld('exit_ufrc_Cu', exit_ufrc, mix, lchnk) 
     call outfld('exit_wtw_Cu', exit_wtw, mix, lchnk) 
     call outfld('exit_drycore_Cu', exit_drycore, mix, lchnk) 
     call outfld('exit_wu_Cu', exit_wu, mix, lchnk) 
     call outfld('exit_cufilter_Cu', exit_cufilter, mix, lchnk) 
     call outfld('exit_kinv1_Cu', exit_kinv1, mix, lchnk) 
     call outfld('exit_rei_Cu', exit_rei, mix, lchnk) 

     call outfld('limit_shcu_Cu', limit_shcu, mix, lchnk) 
     call outfld('limit_negcon_Cu', limit_negcon, mix, lchnk) 
     call outfld('limit_ufrc_Cu', limit_ufrc, mix, lchnk) 
     call outfld('limit_ppen_Cu', limit_ppen, mix, lchnk) 
     call outfld('limit_emf_Cu', limit_emf, mix, lchnk) 
     call outfld('limit_cinlcl_Cu', limit_cinlcl, mix, lchnk) 
     call outfld('limit_cin_Cu', limit_cin, mix, lchnk) 
     call outfld('limit_cbmf_Cu', limit_cbmf, mix, lchnk) 
     call outfld('limit_rei_Cu', limit_rei, mix, lchnk) 
     call outfld('ind_delcin_Cu', ind_delcin, mix, lchnk) 

    return

  end subroutine compute_mcshallow

  ! ------------------------------ !
  !                                ! 
  ! Beginning of subroutine blocks !
  !                                !
  ! ------------------------------ !

  subroutine getbuoy(pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin)
  ! ----------------------------------------------------------- !
  ! Subroutine to calculate integrated CIN [ J/kg = m2/s2 ] and !
  ! 'cinlcl, plfc' if any. Assume 'thv' is linear in each layer !
  ! both for cumulus and environment. Note that this subroutine !
  ! only include positive CIN in calculation - if there are any !
  ! negative CIN, it is assumed to be zero.    This is slightly !
  ! different from 'single_cin' below, where both positive  and !
  ! negative CIN are included.                                  !
  ! ----------------------------------------------------------- !
    real(r8) pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin,frc

    if( thvubot .gt. thv0bot .and. thvutop .gt. thv0top ) then
       plfc = pbot
       return
    elseif( thvubot .le. thv0bot .and. thvutop .le. thv0top ) then 
       cin  = cin - ((thvubot/thv0bot - 1._r8) + (thvutop/thv0top - 1._r8)) * (pbot - ptop)/          &
                    (pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top*exnf(ptop)))
    elseif( thvubot .gt. thv0bot .and. thvutop .le. thv0top ) then 
       frc  = (thvutop/thv0top - 1._r8)/((thvutop/thv0top - 1._r8) - (thvubot/thv0bot - 1._r8))
       cin  = cin - (thvutop/thv0top - 1._r8)*((ptop + frc*(pbot - ptop)) - ptop)/                 &
                    (pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top*exnf(ptop)))
    else            
       frc  = (thvubot/thv0bot - 1._r8)/((thvubot/thv0bot - 1._r8) - (thvutop/thv0top - 1._r8))
       plfc = pbot - frc*(pbot - ptop)
       cin  = cin - (thvubot/thv0bot - 1._r8)*(pbot - plfc)/                                       & 
                    (pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top * exnf(ptop)))
    endif

    return
  end subroutine getbuoy


  function single_cin(pbot,thv0bot,ptop,thv0top,thvubot,thvutop)
  ! ------------------------------------------------------- !
  ! Function to calculate a single layer CIN by summing all ! 
  ! positive and negative CIN.                              !
  ! ------------------------------------------------------- ! 
    real(r8) :: single_cin
    real(r8)    pbot,thv0bot,ptop,thv0top,thvubot,thvutop 

    single_cin = ((1._r8 - thvubot/thv0bot) + (1._r8 - thvutop/thv0top))*(pbot - ptop)/         &
                 (pbot/(r*thv0bot*exnf(pbot)) + ptop/(r*thv0top*exnf(ptop)))
    return
  end function single_cin   


  subroutine conden(p,thl,qt,th,qv,ql,qi,rvls,id_check,qsat)
  ! --------------------------------------------------------------------- !
  ! Calculate thermodynamic properties from a given set of ( p, thl, qt ) !
  ! --------------------------------------------------------------------- !
    implicit none
    real(r8), intent(in)  :: p
    real(r8), intent(in)  :: thl
    real(r8), intent(in)  :: qt
    real(r8), intent(out) :: th
    real(r8), intent(out) :: qv
    real(r8), intent(out) :: ql
    real(r8), intent(out) :: qi
    real(r8), intent(out) :: rvls
    integer , intent(out) :: id_check
    integer , external    :: qsat
    real(r8)              :: tc,temps,t
    real(r8)              :: leff, nu, qc
    integer               :: iteration
    real(r8)              :: es(1)              ! saturation vapor pressure
    real(r8)              :: qs(1)              ! saturation spec. humidity
    real(r8)              :: gam(1)             ! (L/cp)*dqs/dT
    integer               :: status             ! return status of qsat call

    tc     = thl*exnf(p)
    ! This fraction of ice can be refined using 'cldwat_fice' subroutine later.
    nu     = max(min((268._r8 - tc)/20._r8,1.0_r8),0.0_r8)  ! Fraction of ice in the condensate. Func. of T.
    leff   = (1._r8 - nu)*xlv + nu*xls             ! This is an estimate that hopefully speeds convergence

    ! --------------------------------------------------------------------------- !
    ! Below "temps" and "rvls" are just initial guesses for iteration loop below. !
    ! Note that the output "temps" from the below iteration loop is "temperature" !
    ! NOT "liquid temperature".                                                   !
    ! --------------------------------------------------------------------------- !

    temps  = tc
    status = qsat(temps,p,es(1),qs(1),gam(1), 1)
    rvls   = qs(1)

    if(qs(1).ge.qt) then  

      id_check = 0
      qv = qt
      qc = 0._r8
      ql = 0._r8
      qi = 0._r8
      th = tc/exnf(p)

    else 

      do iteration = 1, 10
         temps  = temps + (( tc - temps )*cp/leff + qt - rvls)/               & 
                          (cp/leff + ep2*leff*rvls/r/temps/temps)
         status = qsat(temps,p,es(1),qs(1),gam(1),1)
         rvls   = qs(1)
      end do
         
      qc = max(qt - qs(1),0._r8)
      qv = qt - qc
      ql = qc*(1._r8 - nu)
      qi = nu*qc
      th = temps/exnf(p)

      if( abs((temps-(leff/cp)*qc)-tc).ge.1._r8) then
        id_check = 1
      else
        id_check = 0
      end if

    end if

    return
  end subroutine conden


  subroutine roots(a,b,c,r1,r2,status)
  ! --------------------------------------------------------- !
  ! Subroutine to solve the second order polynomial equation. !
  ! I should check this subroutine later.                     !
  ! --------------------------------------------------------- !
    real(r8), intent(in)  :: a
    real(r8), intent(in)  :: b
    real(r8), intent(in)  :: c
    real(r8), intent(out) :: r1
    real(r8), intent(out) :: r2
    integer , intent(out) :: status
    real(r8)              :: q

    status = 0
    if(a .eq. 0) then                        ! form b*x + c = 0
       if(b .eq. 0) then                     ! failure: c = 0
          status = 1
       else                                  ! b*x + c = 0
          r1 = -c/b
       endif
       r2 = r1
    else
       if(b .eq. 0._r8) then                    ! form a*x**2 + c = 0
          if(a*c .gt. 0._r8) then               ! failure: x**2 = -c/a < 0
             status = 2  
          else                               ! x**2 = -c/a 
             r1 = sqrt(-c/a)
          endif
          r2 = -r1
       else                                  ! form a*x**2 + b*x + c = 0
          if((b**2 - 4._r8*a*c) .lt. 0._r8) then   ! failure, no real roots
             status = 3
          else
             q  = -0.5_r8*(b + sign(1.0_r8,b)*sqrt(b**2 - 4._r8*a*c))
             r1 =  q/a
             r2 =  c/q
          endif
       endif
    endif

    return
  end subroutine roots

  
  function slope(mkx,field,p0)
  ! ------------------------------------------------------------------ !
  ! Function performing profile reconstruction of conservative scalars !
  ! in each layer. This is identical to profile reconstruction used in !
  ! UW-PBL scheme but from bottom to top layer here.     At the lowest !
  ! layer near to surface, slope is defined using the two lowest layer !
  ! mid-point values. I checked this subroutine and it is correct.     !
  ! ------------------------------------------------------------------ !
    real(r8)             :: slope(mkx)
    integer,  intent(in) :: mkx
    real(r8), intent(in) :: field(mkx)
    real(r8), intent(in) :: p0(mkx)
    
    real(r8)             :: below
    real(r8)             :: above
    integer              :: k

    below = (field(2) - field(1))/(p0(2) - p0(1))
    do k = 2, mkx
       above = (field(k) - field(k-1))/(p0(k) - p0(k-1))
       if (above .gt. 0._r8) then
          slope(k-1) = max(0._r8,min(above,below))
       else 
          slope(k-1) = min(0._r8,max(above,below))
       end if
       below = above
    end do
    slope(mkx) = slope(mkx-1)

    return
  end function slope


  function qsinvert(qt,thl,psfc,qsat)
  ! ----------------------------------------------------------------- !
  ! Function calculating saturation pressure ps (or pLCL) from qt and !
  ! thl ( liquid potential temperature,  NOT liquid virtual potential ! 
  ! temperature) by inverting Bolton formula. I should check later if !
  ! current use of 'leff' instead of 'xlv' here is reasonable or not. !
  ! ----------------------------------------------------------------- !
    real(r8)          :: qsinvert    
    real(r8)             qt, thl, psfc
    real(r8)             ps, Pis, Ts, err, dlnqsdT, dTdPis
    real(r8)             dPisdps, dlnqsdps, derrdps, dps 
    real(r8)             Ti, rhi, TLCL, PiLCL, psmin, dpsmax
    integer              i
    integer, external :: qsat
    real(r8)          :: es(1)                     ! saturation vapor pressure
    real(r8)          :: qs(1)                     ! saturation spec. humidity
    real(r8)          :: gam(1)                    ! (L/cp)*dqs/dT
    integer           :: status                    ! return status of qsat call
    real(r8)          :: leff, nu

    psmin = 100._r8*100 ! Default saturation pressure [Pa] if iteration does not converge
    dpsmax = 1._r8      ! Tolerance [Pa] for convergence of iteration

    ! ------------------------------------ !
    ! Calculate best initial guess of pLCL !
    ! ------------------------------------ !

      Ti       =  thl*(psfc/p00)**rovcp
      status   =  qsat(Ti,psfc,es(1),qs(1),gam(1),1)
      rhi      =  qt/qs(1)      
      if( rhi.le.0.01_r8 ) then
        write(iulog,*) 'Source air is too dry and pLCL is set to psmin in mcshallow.F90' 
        qsinvert = psmin
        return
      end if
      TLCL     =  55._r8 + 1._r8/(1._r8/(Ti-55._r8)-log(rhi)/2840._r8); ! Bolton's formula. MWR.1980.Eq.(22)
      PiLCL    =  TLCL/thl
      ps       =  p00*(PiLCL)**(1._r8/rovcp)

    do i = 1, 10
      Pis      =  (ps/p00)**rovcp
      Ts       =  thl*Pis
      status   =  qsat(Ts,ps,es(1),qs(1),gam(1),1)
      err      =  qt - qs(1)
      nu       =  max(min((268._r8 - Ts)/20._r8,1.0_r8),0.0_r8)        
      leff     =  (1._r8 - nu)*xlv + nu*xls                   
      dlnqsdT  =  gam(1)*(cp/leff)/qs(1)
      dTdPis   =  thl
      dPisdps  =  rovcp*Pis/ps 
      dlnqsdps = -1._r8/(ps - (1._r8 - ep2)*es(1))
      derrdps  = -qs(1)*(dlnqsdT * dTdPis * dPisdps + dlnqsdps)
      dps      = -err/derrdps
      ps       =  ps + dps
      if(ps .lt. 0._r8 ) then
         write(iulog,*) 'pLCL iteration is negative and set to psmin in mcshallow.F90' 
         qsinvert = psmin
         return    
      end if
      if( abs(dps) .le. dpsmax ) then
         qsinvert = ps
         return
      end if
    end do
    write(iulog,*) 'pLCL does not converge and is set to psmin in mcshallow.F90' 
    qsinvert = psmin
    return
  end function qsinvert


  real(r8) function compute_alpha(del_CIN,ke)
  ! ------------------------------------------------ !
  ! Subroutine to compute proportionality factor for !
  ! implicit CIN calculation.                        !   
  ! ------------------------------------------------ !
    real(r8) :: del_CIN, ke
    real(r8) :: x0, x1
    integer  :: iteration

    x0 = 0._r8
    do iteration = 1, 10
       x1 = x0 - (exp(-x0*ke*del_CIN) - x0)/(-ke*del_CIN*exp(-x0*ke*del_CIN) - 1._r8)
       x0 = x1
    end do
    compute_alpha = x0

    return

  end function compute_alpha


  real(r8) function compute_mumin2(mulcl,rmaxfrac,mulow)
  ! --------------------------------------------------------- !
  ! Subroutine to compute critical 'mu' (normalized CIN) such ! 
  ! that updraft fraction at the LCL is equal to 'rmaxfrac'.  !
  ! --------------------------------------------------------- !  
    real(r8) :: mulcl, rmaxfrac, mulow
    real(r8) :: x0, x1, ex, ef, exf, f, fs
    integer  :: iteration

    x0 = mulow
    do iteration = 1, 10
       ex = exp(-x0**2)
       ef = erfc(x0)
       ! if(x0.ge.3._r8) then
       !    compute_mumin2 = 3._r8 
       !    goto 20
       ! endif 
       exf = ex/ef
       f  = 0.5_r8*exf**2 - 0.5_r8*(ex/2/rmaxfrac)**2 - (mulcl*2.5066_r8/2)**2
       fs = (2*exf**2)*(exf/sqrt(3.141592_r8)-x0) + (0.5_r8*x0*ex**2)/(rmaxfrac**2)
       x1 = x0 - f/fs     
       x0 = x1
    end do
    compute_mumin2 = x0

 20 return

  end function compute_mumin2


  real(r8) function compute_ppen(wtwb,D,bogbot,bogtop,rho0j,dpen)
  ! ----------------------------------------------------------- !
  ! Subroutine to compute critical 'ppen[Pa]<0' ( pressure dis. !
  ! from 'ps0(kpen-1)' to the cumulus top where cumulus updraft !
  ! vertical velocity is exactly zero ) by considering exact    !
  ! non-zero fer(kpen).                                         !  
  ! ----------------------------------------------------------- !  
    real(r8) :: wtwb, D, bogbot, bogtop, rho0j, dpen
    real(r8) :: x0, x1, f, fs, SB, s00
    integer  :: iteration

    ! Buoyancy slope
      SB = (bogtop-bogbot)/dpen
    ! Sign of slope, 'f' at x = 0
    ! If 's00>0', 'w' increases with height.
      s00 = bogbot/rho0j - D*wtwb

    if( D*dpen.lt.1.e-8 ) then
      if( s00.ge.0._r8 ) then
        x0 = dpen       
      else
        x0 = max(0._r8,min(dpen,-0.5_r8*wtwb/s00))
      endif
    else
      if( s00.ge.0._r8 ) then
          x0 = dpen
      else 
          x0 = 0._r8
      endif
      do iteration = 1, 5
         f  = exp(-2._r8*D*x0)*(wtwb-(bogbot-SB/(2._r8*D))/(D*rho0j)) + &
                            (SB*x0+bogbot-SB/(2._r8*D))/(D*rho0j)
         fs = -2._r8*D*exp(-2._r8*D*x0)*(wtwb-(bogbot-SB/(2._r8*D))/(D*rho0j)) + &
                            (SB)/(D*rho0j)
         x1 = x0 - f/fs     
         x0 = x1
      end do

    endif    

    compute_ppen = -max(0._r8,min(dpen,x0))

  end function compute_ppen


  subroutine fluxbelowinv(cbmf,ps0,mkx,kinv,dt,xsrc,xmean,xtopin,xbotin,xflx)   
  ! ------------------------------------------------------------------------- !
  ! Subroutine to calculate turbulent fluxes at and below 'kinv-1' interfaces.!
  ! Check in the main program such that input 'cbmf' should not be zero.      !  
  ! If the reconstructed inversion height does not go down below the 'kinv-1' !
  ! interface, then turbulent flux at 'kinv-1' interface  is simply a product !
  ! of 'cmbf' and 'qtsrc-xbot' where 'xbot' is the value at the top interface !
  ! of 'kinv-1' layer. This flux is linearly interpolated down to the surface !
  ! assuming turbulent fluxes at surface are zero. If reconstructed inversion !
  ! height goes down below the 'kinv-1' interface, subsidence warming &drying !
  ! measured by 'xtop-xbot', where  'xtop' is the value at the base interface !
  ! of 'kinv+1' layer, is added ONLY to the 'kinv-1' layer, using appropriate !
  ! mass weighting ( rpinv and rcbmf, or rr = rpinv / rcbmf ) between current !
  ! and next provisional time step. Also impose a limiter to enforce outliers !
  ! of thermodynamic variables in 'kinv' layer  to come back to normal values !
  ! at the next step.                                                         !
  ! ------------------------------------------------------------------------- !            
    integer,  intent(in)                     :: mkx, kinv 
    real(r8), intent(in)                     :: cbmf, dt, xsrc, xmean, xtopin, xbotin
    real(r8), intent(in),  dimension(0:mkx)  :: ps0
    real(r8), intent(out), dimension(0:mkx)  :: xflx  
    integer k
    real(r8) rcbmf, rpeff, dp, rr, pinv_eff, xtop, xbot, pinv, xtop_ori, xbot_ori

    xflx(0:mkx) = 0._r8
    dp = ps0(kinv-1)-ps0(kinv) 
    
    if( abs(xbotin-xtopin).le.1.e-33_r8 ) then
      xbot = xbotin - 1.e-20_r8
      xtop = xtopin + 1.e-20_r8
    else
      xbot = xbotin
      xtop = xtopin      
    endif
    ! -------------------------------------- !
    ! Compute reconstructed inversion height !
    ! -------------------------------------- !
    xtop_ori = xtop
    xbot_ori = xbot
    rcbmf = ( cbmf * g * dt ) / dp ! Can be larger than 1 : 'OK'      
    rpeff = ( xmean - xtop ) / ( xbot - xtop ) 
    rpeff = min( max(0._r8,rpeff), 1._r8 )  ! As of this, 0<= rpeff <= 1   
    if(rpeff.eq.0.or.rpeff.eq.1) then
       xbot = xmean
       xtop = xmean
    endif
    ! Below two commented-out lines are the old code replacing the above 'if' block.   
    ! if(rpeff.eq.1) xbot = xmean
    ! if(rpeff.eq.0) xtop = xmean    
    rr = rpeff / rcbmf
    pinv = ps0(kinv-1) - rpeff * dp ! "pinv" before detraining mass
    pinv_eff = ps0(kinv-1) + ( rcbmf - rpeff ) * dp ! Effective "pinv" after detraining mass
    ! -------------------------------------------------------------------- !
    ! Compute turbulent fluxes.                                            !
    ! Below two cases exactly converges at 'kinv-1' interface when rr = 1._r8 !
    ! -------------------------------------------------------------------- !
    do k = 0, kinv - 1
       xflx(k) = cbmf * ( xsrc - xbot ) * ( ps0(0) - ps0(k) ) / ( ps0(0) - pinv )
    end do
    if( rr.le.1 ) then
       xflx(kinv-1) =  xflx(kinv-1) - ( 1._r8 - rr ) * cbmf * ( xtop_ori - xbot_ori )
    endif

    return
  end subroutine fluxbelowinv

  ! ------------------------ !
  !                          ! 
  ! End of subroutine blocks !
  !                          !
  ! ------------------------ !


  end module mcshallow
