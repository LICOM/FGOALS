
module radsw
!----------------------------------------------------------------------- 
! 
! Purpose: Solar radiation calculations.
!
!-----------------------------------------------------------------------
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use abortutils,      only: endrun
use cam_history,     only: outfld
use scamMod,         only: single_column,scm_crm_mode,have_asdir, &
                           asdirobs, have_asdif, asdifobs, have_aldir, &
                           aldirobs, have_aldif, aldifobs
use cam_logfile,     only: iulog
use parrrsw,         only: nbndsw, ngptsw
use rrtmg_sw_init,   only: rrtmg_sw_ini
use rrtmg_sw_rad,    only: rrtmg_sw
use perf_mod,        only: t_startf, t_stopf
use radconstants,    only: idx_sw_diag

implicit none

private
save

real(r8), public :: wavmin(nbndsw) = &  ! Min wavelength (micro-meters) of interval
           (/ 3.077_r8, 2.500_r8, 2.150_r8, 1.942_r8, 1.626_r8, 1.299_r8, 1.242_r8, &
              0.778_r8, 0.625_r8, 0.442_r8, 0.345_r8, 0.263_r8, 0.200_r8, 3.846_r8/)

real(r8), public :: wavmax(nbndsw) = &  ! Max wavelength (micro-meters) of interval
           (/ 3.846_r8, 3.077_r8, 2.500_r8, 2.150_r8, 1.942_r8, 1.626_r8, 1.299_r8, &
              1.242_r8, 0.778_r8, 0.625_r8, 0.442_r8, 0.345_r8, 0.263_r8,12.195_r8/)

real(r8) :: fractional_solar_irradiance(1:nbndsw) ! fraction of solar irradiance in each band
real(r8) :: solar_band_irrad(1:nbndsw) ! rrtmg-assumed solar irradiance in each sw band

! Public methods

public ::&
   radsw_init,      &! initialize constants
   rad_rrtmg_sw      ! driver for solar radiation code

!===============================================================================
CONTAINS
!===============================================================================

subroutine rad_rrtmg_sw(lchnk   ,ncol    ,                        &
                    E_pint    ,E_pmid    ,E_tint    ,E_tnm      , &
                    E_h2o     ,E_co2mmr  ,E_o3mmr ,E_ch4mmr ,E_o2mmr, E_n2ommr, &
                    E_cld     , &
                    E_aer_tau, E_aer_tau_w, E_aer_tau_w_g, E_aer_tau_w_f,   &
                    eccf      ,E_coszrs  ,solin   , sfac, &
                    E_asdir   ,E_asdif   ,E_aldir   ,E_aldif    , &
                    qrs       ,qrsc      ,fsnt      ,fsntc   ,fsntoa, fsutoa, &
                    fsntoac ,fsnirtoa,fsnrtoac,fsnrtoaq,fsns    , &
                    fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                    solsd   ,solld   , fns     ,fcns            , &
                    Nday    ,Nnite   ,IdxDay  ,IdxNite          ,&
                    E_cld_tau, E_cld_tau_w, E_cld_tau_w_g, E_cld_tau_w_f,  &
                    old_convert, ancientmethod)


!-----------------------------------------------------------------------
! 
! Purpose: 
! Solar radiation code
! 
! Method: 
! mji/rrtmg
! RRTMG, two-stream, with McICA
! 
! Divides solar spectrum into 14 intervals from 0.2-12.2 micro-meters.
! solar flux fractions specified for each interval. allows for
! seasonally and diurnally varying solar input.  Includes molecular,
! cloud, aerosol, and surface scattering, along with h2o,o3,co2,o2,cloud, 
! and surface absorption. Computes delta-eddington reflections and
! transmissions assuming homogeneously mixed layers. Adds the layers 
! assuming scattering between layers to be isotropic, and distinguishes 
! direct solar beam from scattered radiation.
! 
! Longitude loops are broken into 1 or 2 sections, so that only daylight
! (i.e. coszrs > 0) computations are done.
! 
! Note that an extra layer above the model top layer is added.
! 
! mks units are used.
! 
! Special diagnostic calculation of the clear sky surface and total column
! absorbed flux is also done for cloud forcing diagnostics.
! 
!-----------------------------------------------------------------------

   use cmparray_mod,        only: CmpDayNite, ExpDayNite
   use phys_control,        only: phys_getopts
   use mcica_subcol_gen_sw, only: mcica_subcol_sw
   use physconst,           only: cpair


   ! Minimum cloud amount (as a fraction of the grid-box area) to 
   ! distinguish from clear sky
   real(r8), parameter :: cldmin = 1.0e-80_r8

   ! Decimal precision of cloud amount (0 -> preserve full resolution;
   ! 10^-n -> preserve n digits of cloud amount)
   real(r8), parameter :: cldeps = 0.0_r8

   ! Input arguments
   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: ncol              ! number of atmospheric columns

   integer, intent(in) :: Nday                      ! Number of daylight columns
   integer, intent(in) :: Nnite                     ! Number of night columns
   integer, intent(in), dimension(pcols) :: IdxDay  ! Indicies of daylight coumns
   integer, intent(in), dimension(pcols) :: IdxNite ! Indicies of night coumns


   real(r8), intent(in) :: E_pmid(pcols,pver)  ! Level pressure (Pascals)
   real(r8), intent(in) :: E_pint(pcols,pverp) ! Interface pressure (Pascals)

   real(r8), intent(in) :: E_tnm(pcols,pver)   ! Level temperature
   real(r8), intent(in) :: E_tint(pcols,pverp) ! Interface temperature
   real(r8), intent(in) :: E_h2o(pcols,pver)   ! Specific humidity

   real(r8), intent(in) :: E_o3mmr(pcols,pver) ! Ozone mass mixing ratio

   real(r8), intent(in) :: E_cld(pcols,pver)    ! Fractional cloud cover

   real(r8), intent(in) :: E_aer_tau    (pcols, 0:pver, nbndsw)      ! aerosol optical depth
   real(r8), intent(in) :: E_aer_tau_w  (pcols, 0:pver, nbndsw)      ! aerosol OD * ssa
   real(r8), intent(in) :: E_aer_tau_w_g(pcols, 0:pver, nbndsw)      ! aerosol OD * ssa * asm
   real(r8), intent(in) :: E_aer_tau_w_f(pcols, 0:pver, nbndsw)      ! aerosol OD * ssa * fwd

   real(r8), intent(in) :: eccf               ! Eccentricity factor (1./earth-sun dist^2)
   real(r8), intent(in) :: E_coszrs(pcols)    ! Cosine solar zenith angle
   real(r8), intent(in) :: E_asdir(pcols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: E_aldir(pcols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: E_asdif(pcols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
   real(r8), intent(in) :: E_aldif(pcols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad

   real(r8), intent(in) :: E_co2mmr(pcols,pver)    ! co2 mass mixing ratio
   real(r8), intent(in) :: E_ch4mmr(pcols,pver)    ! ch4 mass mixing ratio
   real(r8), intent(in) :: E_o2mmr(pcols,pver)     ! o2  mass mixing ratio
   real(r8), intent(in) :: E_n2ommr(pcols,pver)    ! n2o column mean mmr
   real(r8), intent(in) :: sfac(nbndsw)            ! factor to account for solar variability in each band 

   real(r8), optional, intent(in) :: E_cld_tau    (nbndsw, pcols, pver)      ! cloud optical depth
   real(r8), optional, intent(in) :: E_cld_tau_w  (nbndsw, pcols, pver)      ! cloud optical 
   real(r8), optional, intent(in) :: E_cld_tau_w_g(nbndsw, pcols, pver)      ! cloud optical 
   real(r8), optional, intent(in) :: E_cld_tau_w_f(nbndsw, pcols, pver)      ! cloud optical 
   logical, optional, intent(in) :: old_convert
   logical, optional, intent(in) :: ancientmethod

   ! Output arguments

   real(r8), intent(out) :: solin(pcols)     ! Incident solar flux
   real(r8), intent(out) :: qrs (pcols,pver) ! Solar heating rate
   real(r8), intent(out) :: qrsc(pcols,pver) ! Clearsky solar heating rate
   real(r8), intent(out) :: fsns(pcols)      ! Surface absorbed solar flux
   real(r8), intent(out) :: fsnt(pcols)      ! Total column absorbed solar flux
   real(r8), intent(out) :: fsntoa(pcols)    ! Net solar flux at TOA
   real(r8), intent(out) :: fsutoa(pcols)    ! Upward solar flux at TOA
   real(r8), intent(out) :: fsds(pcols)      ! Flux shortwave downwelling surface

   real(r8), intent(out) :: fsnsc(pcols)     ! Clear sky surface absorbed solar flux
   real(r8), intent(out) :: fsdsc(pcols)     ! Clear sky surface downwelling solar flux
   real(r8), intent(out) :: fsntc(pcols)     ! Clear sky total column absorbed solar flx
   real(r8), intent(out) :: fsntoac(pcols)   ! Clear sky net solar flx at TOA
   real(r8), intent(out) :: sols(pcols)      ! Direct solar rad on surface (< 0.7)
   real(r8), intent(out) :: soll(pcols)      ! Direct solar rad on surface (>= 0.7)
   real(r8), intent(out) :: solsd(pcols)     ! Diffuse solar rad on surface (< 0.7)
   real(r8), intent(out) :: solld(pcols)     ! Diffuse solar rad on surface (>= 0.7)
   real(r8), intent(out) :: fsnirtoa(pcols)  ! Near-IR flux absorbed at toa
   real(r8), intent(out) :: fsnrtoac(pcols)  ! Clear sky near-IR flux absorbed at toa
   real(r8), intent(out) :: fsnrtoaq(pcols)  ! Net near-IR flux at toa >= 0.7 microns

   real(r8), intent(out) :: fns(pcols,pverp)   ! net flux at interfaces
   real(r8), intent(out) :: fcns(pcols,pverp)  ! net clear-sky flux at interfaces

   !---------------------------Local variables-----------------------------

   ! Local and reordered copies of the intent(in) variables

   real(r8) :: pmid(pcols,pver)    ! Level pressure (Pascals)
   real(r8) :: pint(pcols,pverp)   ! Interface pressure (Pascals)
   real(r8) :: pmidmb(pcols,pver)  ! Level pressure (hPa)
   real(r8) :: pintmb(pcols,pverp) ! Interface pressure (hPa)

   real(r8) :: tnm(pcols,pver)   ! Level temperature
   real(r8) :: tint(pcols,pverp) ! Interface temperature
   real(r8) :: h2o(pcols,pver)   ! Specific humidity

   real(r8) :: o3mmr(pcols,pver) ! Ozone mass mixing ratio

   real(r8) :: cld(pcols,pver)    ! Fractional cloud cover
   real(r8) :: cicewp(pcols,pver) ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,pver) ! in-cloud cloud liquid water path
   real(r8) :: rel(pcols,pver)    ! Liquid effective drop size (microns)
   real(r8) :: rei(pcols,pver)    ! Ice effective drop size (microns)

   real(r8) :: coszrs(pcols)    ! Cosine solar zenith angle
   real(r8) :: asdir(pcols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
   real(r8) :: aldir(pcols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
   real(r8) :: asdif(pcols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
   real(r8) :: aldif(pcols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad

   real(r8) :: co2mmr(pcols,pver)   ! co2 mass mixing ratio
   real(r8) :: ch4mmr(pcols,pver)   ! ch4 mass mixing ratio
   real(r8) :: o2mmr(pcols,pver)    ! o2  mass mixing ratio
   real(r8) :: n2ommr(pcols,pver)   ! n2o mass mixing ratio

   real(r8) :: h2ovmr(pcols,pver)   ! h2o volume mixing ratio
   real(r8) :: o3vmr(pcols,pver)    ! o3 volume mixing ratio
   real(r8) :: co2vmr(pcols,pver)   ! co2 volume mixing ratio 
   real(r8) :: ch4vmr(pcols,pver)   ! ch4 volume mixing ratio 
   real(r8) :: o2vmr(pcols,pver)    ! o2  volume mixing ratio 
   real(r8) :: n2ovmr(pcols,pver)   ! n2o volume mixing ratio 
   real(r8) :: tsfc(pcols)          ! surface temperature

   ! Set molecular weight ratios (for converting mmr to vmr)
   !  e.g. h2ovmr = h2ommr * amdw)
   real(r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
   real(r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
   real(r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
   real(r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
   real(r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
   real(r8), parameter :: amdo2 = 0.905140_r8   ! Molecular weight of dry air / oxygen

   integer :: inflgsw               ! flag for cloud parameterization method
   integer :: iceflgsw              ! flag for ice cloud parameterization method
   integer :: liqflgsw              ! flag for liquid cloud parameterization method
   integer :: icld                  ! Flag for cloud overlap method
                                    ! 0=clear, 1=random, 2=maximum/random, 3=maximum
   integer :: dyofyr                ! Set to day of year for Earth/Sun distance calculation in
                                    ! rrtmg_sw, or pass in adjustment directly into adjes
   real(r8) :: solvar(nbndsw)       ! solar irradiance variability in each band

   integer, parameter :: nsubcsw = ngptsw           ! rrtmg_sw g-point (quadrature point) dimension
   integer :: permuteseed                           ! permute seed for sub-column generator

   real(r8) :: diagnostic_od(pcols, pver)           ! cloud optical depth - diagnostic temp variable

   real(r8) :: tauc_sw(nbndsw, pcols, pver)         ! cloud optical depth
   real(r8) :: ssac_sw(nbndsw, pcols, pver)         ! cloud single scat. albedo
   real(r8) :: asmc_sw(nbndsw, pcols, pver)         ! cloud asymmetry parameter
   real(r8) :: fsfc_sw(nbndsw, pcols, pver)         ! cloud forward scattering fraction

   real(r8) :: tau_aer_sw(pcols, pver, nbndsw)      ! aer optical depth
   real(r8) :: ssa_aer_sw(pcols, pver, nbndsw)      ! aer single scat. albedo
   real(r8) :: asm_aer_sw(pcols, pver, nbndsw)      ! aer asymmetry parameter

   real(r8) :: cld_stosw(nsubcsw, pcols, pver)      ! stochastic cloud fraction
   real(r8) :: rei_stosw(pcols, pver)               ! stochastic ice particle size 
   real(r8) :: rel_stosw(pcols, pver)               ! stochastic liquid particle size
   real(r8) :: cicewp_stosw(nsubcsw, pcols, pver)   ! stochastic cloud ice water path
   real(r8) :: cliqwp_stosw(nsubcsw, pcols, pver)   ! stochastic cloud liquid wter path
   real(r8) :: tauc_stosw(nsubcsw, pcols, pver)     ! stochastic cloud optical depth (optional)
   real(r8) :: ssac_stosw(nsubcsw, pcols, pver)     ! stochastic cloud single scat. albedo (optional)
   real(r8) :: asmc_stosw(nsubcsw, pcols, pver)     ! stochastic cloud asymmetry parameter (optional)
   real(r8) :: fsfc_stosw(nsubcsw, pcols, pver)     ! stochastic cloud forward scattering fraction (optional)

   real(r8), parameter :: dps = 1._r8/86400._r8 ! Inverse of seconds per day
 
   real(r8) :: swuflx(pcols,pverp+1)       ! Total sky shortwave upward flux (W/m2)
   real(r8) :: swdflx(pcols,pverp+1)       ! Total sky shortwave downward flux (W/m2)
   real(r8) :: swhr(pcols,pver+1)          ! Total sky shortwave radiative heating rate (K/d)
   real(r8) :: swuflxc(pcols,pverp+1)      ! Clear sky shortwave upward flux (W/m2)
   real(r8) :: swdflxc(pcols,pverp+1)      ! Clear sky shortwave downward flux (W/m2)
   real(r8) :: swhrc(pcols,pver+1)         ! Clear sky shortwave radiative heating rate (K/d)

   real(r8) :: dirdnuv(pcols,pverp+1)       ! Direct downward shortwave flux, UV/vis
   real(r8) :: difdnuv(pcols,pverp+1)       ! Diffuse downward shortwave flux, UV/vis
   real(r8) :: dirdnir(pcols,pverp+1)       ! Direct downward shortwave flux, near-IR
   real(r8) :: difdnir(pcols,pverp+1)       ! Diffuse downward shortwave flux, near-IR

   ! Added for net near-IR diagnostic
   real(r8) :: ninflx(pcols,pverp+1)        ! Net shortwave flux, near-IR
   real(r8) :: ninflxc(pcols,pverp+1)       ! Net clear sky shortwave flux, near-IR

   ! Other

   integer :: i, k, ns       ! indices

   ! Cloud radiative property arrays
   real(r8) :: tauxcl(pcols,0:pver) ! water cloud extinction optical depth
   real(r8) :: tauxci(pcols,0:pver) ! ice cloud extinction optical depth
   real(r8) :: wcl(pcols,0:pver) ! liquid cloud single scattering albedo
   real(r8) :: gcl(pcols,0:pver) ! liquid cloud asymmetry parameter
   real(r8) :: fcl(pcols,0:pver) ! liquid cloud forward scattered fraction
   real(r8) :: wci(pcols,0:pver) ! ice cloud single scattering albedo
   real(r8) :: gci(pcols,0:pver) ! ice cloud asymmetry parameter
   real(r8) :: fci(pcols,0:pver) ! ice cloud forward scattered fraction

   ! Aerosol radiative property arrays
   real(r8) :: tauxar(pcols,0:pver) ! aerosol extinction optical depth
   real(r8) :: wa(pcols,0:pver) ! aerosol single scattering albedo
   real(r8) :: ga(pcols,0:pver) ! aerosol assymetry parameter
   real(r8) :: fa(pcols,0:pver) ! aerosol forward scattered fraction

   ! CRM
   real(r8) :: fus(pcols,pverp)   ! Upward flux (added for CRM)
   real(r8) :: fds(pcols,pverp)   ! Downward flux (added for CRM)
   real(r8) :: fusc(pcols,pverp)  ! Upward clear-sky flux (added for CRM)
   real(r8) :: fdsc(pcols,pverp)  ! Downward clear-sky flux (added for CRM)


   !-----------------------------------------------------------------------
   ! START OF CALCULATION
   !-----------------------------------------------------------------------

   ! Initialize output fields:

   fsds(1:ncol)     = 0.0_r8

   fsnirtoa(1:ncol) = 0.0_r8
   fsnrtoac(1:ncol) = 0.0_r8
   fsnrtoaq(1:ncol) = 0.0_r8

   fsns(1:ncol)     = 0.0_r8
   fsnsc(1:ncol)    = 0.0_r8
   fsdsc(1:ncol)    = 0.0_r8

   fsnt(1:ncol)     = 0.0_r8
   fsntc(1:ncol)    = 0.0_r8
   fsntoa(1:ncol)   = 0.0_r8
   fsutoa(1:ncol)   = 0.0_r8
   fsntoac(1:ncol)  = 0.0_r8

   solin(1:ncol)    = 0.0_r8

   sols(1:ncol)     = 0.0_r8
   soll(1:ncol)     = 0.0_r8
   solsd(1:ncol)    = 0.0_r8
   solld(1:ncol)    = 0.0_r8

   qrs (1:ncol,1:pver) = 0.0_r8
   qrsc(1:ncol,1:pver) = 0.0_r8
   fns(1:ncol,1:pverp) = 0.0_r8
   fcns(1:ncol,1:pverp) = 0.0_r8
   if (single_column.and.scm_crm_mode) then 
      fus(1:ncol,1:pverp) = 0.0_r8
      fds(1:ncol,1:pverp) = 0.0_r8
      fusc(:ncol,:pverp) = 0.0_r8
      fdsc(:ncol,:pverp) = 0.0_r8
   endif

   ! If night everywhere, return:
   if ( Nday == 0 ) then
     return
   endif

   ! Rearrange input arrays
   call CmpDayNite(E_pmid, pmid,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_pint, pint,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call CmpDayNite(E_h2o, h2o,  	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_o3mmr, o3mmr,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_cld, cld,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_co2mmr, co2mmr,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_coszrs, coszrs,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_asdir, asdir,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_aldir, aldir,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_asdif, asdif,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_aldif, aldif,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)

   call CmpDayNite(E_tnm, tnm,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_tint, tint,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call CmpDayNite(E_ch4mmr, ch4mmr,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_o2mmr, o2mmr,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_n2ommr, n2ommr,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)

   ! These fields are no longer input by CAM.
   cicewp = 0.0_r8
   cliqwp = 0.0_r8
   rel = 0.0_r8
   rei = 0.0_r8

   ! Aerosol daylight map
   ! Also convert to optical properties of rrtmg interface, even though
   !   these quantities are later multiplied back together inside rrtmg !
   ! Why does rrtmg use the factored quantities?
   ! There are several different ways this factoring could be done.
   ! Other ways might allow for better optimization
   do ns = 1, nbndsw
      do k  = 1, pver
         do i  = 1, Nday
            if(E_aer_tau_w(IdxDay(i),k,ns) > 1.e-80_r8) then
               asm_aer_sw(i,k,ns) = E_aer_tau_w_g(IdxDay(i),k,ns)/E_aer_tau_w(IdxDay(i),k,ns)
            else
               asm_aer_sw(i,k,ns) = 0._r8
            endif
            if(E_aer_tau(IdxDay(i),k,ns) > 0._r8) then
               ssa_aer_sw(i,k,ns) = E_aer_tau_w(IdxDay(i),k,ns)/E_aer_tau(IdxDay(i),k,ns)
               tau_aer_sw(i,k,ns) = E_aer_tau(IdxDay(i),k,ns)
            else
               ssa_aer_sw(i,k,ns) = 1._r8
               tau_aer_sw(i,k,ns) = 0._r8
            endif
         enddo
      enddo
   enddo

   if (scm_crm_mode) then
      ! overwrite albedos for CRM
      if(have_asdir) asdir = asdirobs(1)
      if(have_asdif) asdif = asdifobs(1)
      if(have_aldir) aldir = aldirobs(1)
      if(have_aldif) aldif = aldifobs(1)
   endif

   ! Define solar incident radiation
   do i = 1, Nday
      solin(i)  = sum(sfac(:)*solar_band_irrad(:)) * eccf * coszrs(i)
   end do

   ! Calculate cloud optical properties here if using CAM method, or if using one of the
   ! methods in RRTMG_SW, then pass in cloud physical properties and zero out cloud optical 
   ! properties here

   ! Zero optional cloud optical property input arrays tauc_sw, ssac_sw, asmc_sw, 
   ! if inputting cloud physical properties to RRTMG_SW
   !tauc_sw(:,:,:) = 0.0_r8
   !ssac_sw(:,:,:) = 1.0_r8
   !asmc_sw(:,:,:) = 0.0_r8
   !fsfc_sw(:,:,:) = 0.0_r8
   !
   ! Or, calculate and pass in CAM cloud shortwave optical properties to RRTMG_SW
   !if (present(old_convert)) print *, 'old_convert',old_convert
   !if (present(ancientmethod)) print *, 'ancientmethod',ancientmethod
   if (present(old_convert))then
      if (old_convert)then ! convert without limits
         do i = 1, Nday
         do k = 1, pver
         do ns = 1, nbndsw
           if (E_cld_tau_w(ns,IdxDay(i),k) > 0._r8) then
              fsfc_sw(ns,i,k)=E_cld_tau_w_f(ns,IdxDay(i),k)/E_cld_tau_w(ns,IdxDay(i),k)
              asmc_sw(ns,i,k)=E_cld_tau_w_g(ns,IdxDay(i),k)/E_cld_tau_w(ns,IdxDay(i),k)
           else
              fsfc_sw(ns,i,k) = 0._r8
              asmc_sw(ns,i,k) = 0._r8
           endif
   
           tauc_sw(ns,i,k)=E_cld_tau(ns,IdxDay(i),k)
           if (tauc_sw(ns,i,k) > 0._r8) then
              ssac_sw(ns,i,k)=E_cld_tau_w (ns,IdxDay(i),k)/tauc_sw(ns,i,k)
           else
              tauc_sw(ns,i,k) = 0._r8
              fsfc_sw(ns,i,k) = 0._r8
              asmc_sw(ns,i,k) = 0._r8
              ssac_sw(ns,i,k) = 1._r8
           endif
         enddo
         enddo
         enddo
      else
         ! eventually, when we are done with archaic versions, This set of code will become the default.
         do i = 1, Nday
         do k = 1, pver
         do ns = 1, nbndsw
           if (E_cld_tau_w(ns,IdxDay(i),k) > 0._r8) then
              fsfc_sw(ns,i,k)=E_cld_tau_w_f(ns,IdxDay(i),k)/max(E_cld_tau_w(ns,IdxDay(i),k), 1.e-80_r8)
              asmc_sw(ns,i,k)=E_cld_tau_w_g(ns,IdxDay(i),k)/max(E_cld_tau_w(ns,IdxDay(i),k), 1.e-80_r8)
           else
              fsfc_sw(ns,i,k) = 0._r8
              asmc_sw(ns,i,k) = 0._r8
           endif
   
           tauc_sw(ns,i,k)=E_cld_tau(ns,IdxDay(i),k)
           if (tauc_sw(ns,i,k) > 0._r8) then
              ssac_sw(ns,i,k)=max(E_cld_tau_w (ns,IdxDay(i),k),1.e-80_r8)/max(tauc_sw(ns,i,k),1.e-80_r8)
           else
              tauc_sw(ns,i,k) = 0._r8
              fsfc_sw(ns,i,k) = 0._r8
              asmc_sw(ns,i,k) = 0._r8
              ssac_sw(ns,i,k) = 1._r8
           endif
         enddo
         enddo
         enddo
      endif
   else
      do i = 1, Nday
      do k = 1, pver
      do ns = 1, nbndsw
        if (E_cld_tau_w(ns,IdxDay(i),k) > 0._r8) then
           fsfc_sw(ns,i,k)=E_cld_tau_w_f(ns,IdxDay(i),k)/max(E_cld_tau_w(ns,IdxDay(i),k), 1.e-80_r8)
           asmc_sw(ns,i,k)=E_cld_tau_w_g(ns,IdxDay(i),k)/max(E_cld_tau_w(ns,IdxDay(i),k), 1.e-80_r8)
        else
           fsfc_sw(ns,i,k) = 0._r8
           asmc_sw(ns,i,k) = 0._r8
        endif

        tauc_sw(ns,i,k)=E_cld_tau(ns,IdxDay(i),k)
        if (tauc_sw(ns,i,k) > 0._r8) then
           ssac_sw(ns,i,k)=max(E_cld_tau_w (ns,IdxDay(i),k),1.e-80_r8)/max(tauc_sw(ns,i,k),1.e-80_r8)
        else
           tauc_sw(ns,i,k) = 0._r8
           fsfc_sw(ns,i,k) = 0._r8
           asmc_sw(ns,i,k) = 0._r8
           ssac_sw(ns,i,k) = 1._r8
        endif
      enddo
      enddo
      enddo
   endif

   ! Call mcica sub-column generator for RRTMG_SW

   ! Call sub-column generator for McICA in radiation
   call t_startf('mcica_subcol_sw')

   ! Select cloud overlap approach (1=random, 2=maximum-random, 3=maximum)
   icld = 2
   ! Set permute seed (must be offset between LW and SW by at least 140 to insure 
   ! effective randomization)
   permuteseed = 1

   call mcica_subcol_sw(lchnk, Nday, pver, icld, permuteseed, pmid, &
      cld, cicewp, cliqwp, rei, rel, tauc_sw, ssac_sw, asmc_sw, fsfc_sw, &
      cld_stosw, cicewp_stosw, cliqwp_stosw, rei_stosw, rel_stosw, &
      tauc_stosw, ssac_stosw, asmc_stosw, fsfc_stosw)

   call t_stopf('mcica_subcol_sw')

   call t_startf('rrtmg_sw')

   ! Call RRTMG_SW for all layers for daylight columns

   ! Select parameterization of cloud ice and liquid optical depths
   ! Use CAM shortwave cloud optical properties directly
   inflgsw = 0 
   iceflgsw = 0
   liqflgsw = 0
   ! Use E&C param for ice to mimic CAM3 for now
   !   inflgsw = 2 
   !   iceflgsw = 1
   !   liqflgsw = 1
   ! Use merged Fu and E&C params for ice 
   !   inflgsw = 2 
   !   iceflgsw = 3
   !   liqflgsw = 1

   ! Set day of year for Earth/Sun distance calculation in rrtmg_sw, or
   ! set to zero and pass E/S adjustment (eccf) directly into array adjes
   dyofyr = 0

   ! Convert incoming water amounts from specific humidity to vmr as needed;
   ! Convert other incoming molecular amounts from mmr to vmr as needed;
   ! Convert pressures from Pa to hPa;
   ! Set surface temperature
   do k = 1, pver
      do i = 1, ncol
         h2ovmr(i,k) = (h2o(i,k) / (1._r8 - h2o(i,k))) * amdw
         o3vmr(i,k) = o3mmr(i,k) * amdo
         co2vmr(i,k) = co2mmr(i,k) * amdc
         ch4vmr(i,k) = ch4mmr(i,k) * amdm
         o2vmr(i,k) = o2mmr(i,k) * amdo2
         n2ovmr(i,k) = n2ommr(i,k) * amdn
         pmidmb(i,k) = pmid(i,k) * 1.e-2_r8
         pintmb(i,k) = pint(i,k) * 1.e-2_r8
      enddo
   enddo
   do i = 1, ncol
      tsfc(i) = tint(i,pverp)
      pintmb(i,pverp) = pint(i,pverp) * 1.e-2_r8
   enddo

   solvar(1:nbndsw) = sfac(1:nbndsw)

   call rrtmg_sw(lchnk, Nday, pver, icld,         &
                 pmidmb, pintmb, tnm, tint, tsfc, &
                 h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, &
                 asdir, asdif, aldir, aldif, &
                 coszrs, eccf, dyofyr, solvar, &
                 inflgsw, iceflgsw, liqflgsw, &
                 cld_stosw, tauc_stosw, ssac_stosw, asmc_stosw, fsfc_stosw, &
                 cicewp_stosw, cliqwp_stosw, rei, rel, &
                 tau_aer_sw, ssa_aer_sw, asm_aer_sw, &
                 swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
                 dirdnuv, dirdnir, difdnuv, difdnir, ninflx, ninflxc)

   ! Flux units are in W/m2 on output from rrtmg_sw and contain output for
   ! extra layer above model top with vertical indexing from bottom to top.
   !
   ! Heating units are in J/kg/s on output from rrtmg_sw and contain output 
   ! for extra layer above model top with vertical indexing from bottom to top.  
   !
   ! Reverse vertical indexing to go from top to bottom for CAM output.

   do i=1,Nday

      ! Set the net absorted shortwave flux at TOA (top of extra layer)
      fsntoa(i) = swdflx(i,pverp+1) - swuflx(i,pverp+1)
      fsutoa(i) = swuflx(i,pverp+1)
      fsntoac(i) = swdflxc(i,pverp+1) - swuflxc(i,pverp+1)

      ! Set net near-IR flux at top of the model
      fsnirtoa(i) = ninflx(i,pverp)
      fsnrtoaq(i) = ninflx(i,pverp)
      fsnrtoac(i) = ninflxc(i,pverp)

      ! Set the net absorbed shortwave flux at the model top level
      fsnt(i) = swdflx(i,pverp) - swuflx(i,pverp)
      fsntc(i) = swdflxc(i,pverp) - swuflxc(i,pverp)
      
      ! Set the downwelling flux at the surface 
      fsds(i) = swdflx(i,1)
      fsdsc(i) = swdflxc(i,1)

      ! Set the net shortwave flux at the surface
      fsns(i) = swdflx(i,1) - swuflx(i,1)
      fsnsc(i) = swdflxc(i,1) - swuflxc(i,1)

      ! Set the UV/vis and near-IR direct and dirruse downward shortwave flux at surface
      sols(i) = dirdnuv(i,1)
      soll(i) = dirdnir(i,1)
      solsd(i) = difdnuv(i,1)
      solld(i) = difdnir(i,1)

      ! Set the net, up and down fluxes at model interfaces
      do k = 1, pverp
         fns(i,k) = swdflx(i,pverp+1-k) - swuflx(i,pverp+1-k)
         fcns(i,k) = swdflxc(i,pverp+1-k) - swuflxc(i,pverp+1-k)
         fus(i,k) = swuflx(i,pverp+1-k)
         fusc(i,k) = swuflxc(i,pverp+1-k)
         fds(i,k) = swdflx(i,pverp+1-k)
         fdsc(i,k) = swdflxc(i,pverp+1-k)
      end do

      ! Set solar heating, reverse layering
      ! Pass shortwave heating to CAM arrays and convert from K/d to J/kg/s
      do k = 1, pver
         qrs(i,k) = swhr(i,pver+1-k)*cpair*dps
         qrsc(i,k) = swhrc(i,pver+1-k)*cpair*dps
      end do

   end do

   call t_stopf('rrtmg_sw')


   ! Rearrange output arrays.
   !
   ! intent(out)
   call ExpDayNite(solin,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(qrs,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call ExpDayNite(qrsc,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call ExpDayNite(fns,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call ExpDayNite(fcns,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call ExpDayNite(fsns,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnt,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsntoa,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsutoa,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsds,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnsc,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsdsc,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsntc,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsntoac,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(sols,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(soll,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(solsd,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(solld,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnirtoa,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnrtoac,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnrtoaq,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)

   !  these outfld calls don't work for spmd only outfield in scm mode (nonspmd)
   if (single_column .and. scm_crm_mode) then 
      ! Following outputs added for CRM
      call ExpDayNite(fus,Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
      call ExpDayNite(fds,Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
      call ExpDayNite(fusc,Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
      call ExpDayNite(fdsc,Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
      call outfld('FUS     ',fus * 1.e-3_r8 ,pcols,lchnk)
      call outfld('FDS     ',fds * 1.e-3_r8 ,pcols,lchnk)
      call outfld('FUSC    ',fusc,pcols,lchnk)
      call outfld('FDSC    ',fdsc,pcols,lchnk)
   endif

end subroutine rad_rrtmg_sw

!-------------------------------------------------------------------------------

subroutine radsw_init()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize various constants for radiation scheme.
!
!-----------------------------------------------------------------------
    use radconstants,  only: get_solar_band_fraction_irrad, get_ref_solar_band_irrad

    ! get the reference fractional solar irradiance in each band
    call get_solar_band_fraction_irrad(fractional_solar_irradiance)
    call get_ref_solar_band_irrad( solar_band_irrad )


   ! Initialize rrtmg_sw
   call rrtmg_sw_ini
 
end subroutine radsw_init


!-------------------------------------------------------------------------------

end module radsw
