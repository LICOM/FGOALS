
module radlw
!----------------------------------------------------------------------- 
! 
! Purpose: Longwave radiation calculations.
!
!-----------------------------------------------------------------------
use shr_kind_mod,      only: r8 => shr_kind_r8
use ppgrid,            only: pcols, pver, pverp
use scamMod,           only: single_column, scm_crm_mode
use parrrtm,           only: nbndlw, ngptlw
use rrtmg_lw_init,     only: rrtmg_lw_ini
use rrtmg_lw_rad,      only: rrtmg_lw
use spmd_utils,        only: masterproc
use perf_mod,          only: t_startf, t_stopf
use cam_logfile,       only: iulog
use abortutils,        only: endrun
use radconstants,      only: nlwbands

implicit none

private
save

! Public methods

public ::&
   radlw_init,   &! initialize constants
   rad_rrtmg_lw   ! driver for longwave radiation code


real(r8), public :: wavenumber1_longwave(nlwbands) ! Longwave spectral band limits (cm-1)
data wavenumber1_longwave &
    /10.,350.,500.,630.,700.,820.,980.,1080.,1180.,1390.,1480.,1800.,2080.,2250.,2390.,2600./
real(r8), public :: wavenumber2_longwave(nlwbands) ! Longwave spectral band limits (cm-1)
data wavenumber2_longwave &
    /350.,500.,630.,700.,820.,980.,1080.,1180.,1390.,1480.,1800.,2080., 2250.,2390.,2600.,3250./

! Private data
integer :: ntoplw    ! top level to solve for longwave cooling

!===============================================================================
CONTAINS
!===============================================================================

subroutine rad_rrtmg_lw(lchnk   ,ncol    ,                   &
                        tint    ,tnm     ,qnm     ,co2mmr  ,o3, &
                        pmid    ,pint    ,aer_lw_abs       , &
                        n2o     ,ch4     ,o2      ,cfc11   ,cfc12   , &
                        cld     ,tauc_lw , &
                        qrl     ,qrlc    , &
                        flns    ,flnt    ,flnsc   ,flntc   ,flwds, &
                        flut    ,flutc   ,fnl     ,fcnl, fldsc)

!-----------------------------------------------------------------------
   use cam_history,         only: outfld
   use mcica_subcol_gen_lw, only: mcica_subcol_lw
   use physconst,           only: cpair

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: tint(pcols,pverp)    ! Interface temperature
!
! Input arguments which are only passed to other routines
!
   real(r8), intent(in) :: tnm(pcols,pver)      ! Level temperature
   real(r8), intent(in) :: qnm(pcols,pver)      ! Level moisture field
   real(r8), intent(in) :: o3(pcols,pver)       ! ozone mass mixing ratio
   real(r8), intent(in) :: pmid(pcols,pver)     ! Level pressure (Pascals)
   real(r8), intent(in) :: pint(pcols,pverp)    ! Model interface pressure (Pascals)

   real(r8), intent(in) :: aer_lw_abs (pcols,pver,nbndlw) ! aerosol absorption optics depth (LW)

   real(r8), intent(in) :: co2mmr(pcols,pver)   ! co2 mass mixing ratio
   real(r8), intent(in) :: n2o(pcols,pver)      ! nitrous oxide mass mixing ratio
   real(r8), intent(in) :: ch4(pcols,pver)      ! methane mass mixing ratio
   real(r8), intent(in) :: o2(pcols,pver)       ! o2 mass mixing ratio
   real(r8), intent(in) :: cfc11(pcols,pver)    ! cfc11 mass mixing ratio
   real(r8), intent(in) :: cfc12(pcols,pver)    ! cfc12 mass mixing ratio
   real(r8), intent(in) :: cld(pcols,pver)      ! Cloud cover
   real(r8), intent(in) :: tauc_lw(nbndlw,pcols,pver)   ! Cloud longwave optical depth by band

!
! Output arguments
!
   real(r8), intent(out) :: qrl (pcols,pver)     ! Longwave heating rate
   real(r8), intent(out) :: qrlc(pcols,pver)     ! Clearsky longwave heating rate
   real(r8), intent(out) :: flns(pcols)          ! Surface cooling flux
   real(r8), intent(out) :: flnt(pcols)          ! Net outgoing flux
   real(r8), intent(out) :: flut(pcols)          ! Upward flux at top of model
   real(r8), intent(out) :: flnsc(pcols)         ! Clear sky surface cooing
   real(r8), intent(out) :: flntc(pcols)         ! Net clear sky outgoing flux
   real(r8), intent(out) :: flutc(pcols)         ! Upward clear-sky flux at top of model
   real(r8), intent(out) :: flwds(pcols)         ! Down longwave flux at surface
   real(r8), intent(out) :: fldsc(pcols)         ! Down longwave clear flux at surface
   real(r8), intent(out) :: fcnl(pcols,pverp)    ! clear sky net flux at interfaces
   real(r8), intent(out) :: fnl(pcols,pverp)     ! net flux at interfaces
!
!---------------------------Local variables-----------------------------
!
   integer :: i, k, nbnd         ! indices

   real(r8) ful(pcols,pverp)     ! Total upwards longwave flux
   real(r8) fsul(pcols,pverp)    ! Clear sky upwards longwave flux
   real(r8) fdl(pcols,pverp)     ! Total downwards longwave flux
   real(r8) fsdl(pcols,pverp)    ! Clear sky downwards longwv flux

   integer inflglw               ! Flag for cloud parameterization method
   integer iceflglw              ! Flag for ice cloud param method
   integer liqflglw              ! Flag for liquid cloud param method
   integer icld                  ! Flag for cloud overlap method
                                 ! 0=clear, 1=random, 2=maximum/random, 3=maximum

   real(r8) h2ovmr(pcols,pver)   ! h2o volume mixing ratio
   real(r8) o3vmr(pcols,pver)    ! o3 volume mixing ratio
   real(r8) co2vmr(pcols,pver)   ! co2 volume mixing ratio 
   real(r8) ch4vmr(pcols,pver)   ! ch4 volume mixing ratio 
   real(r8) o2vmr(pcols,pver)    ! o2  volume mixing ratio 
   real(r8) n2ovmr(pcols,pver)   ! n2o volume mixing ratio 
   real(r8) cfc11vmr(pcols,pver) ! cfc11 volume mixing ratio
   real(r8) cfc12vmr(pcols,pver) ! cfc12 volume mixing ratio
   real(r8) cfc22vmr(pcols,pver) ! cfc22 volume mixing ratio
   real(r8) ccl4vmr(pcols,pver)  ! ccl4 volume mixing ratio
   real(r8) tsfc(pcols)          ! surface temperature
   real(r8) emis(pcols,nbndlw)   ! surface emissivity

   real(r8) pmidmb(pcols,pver)     ! Level pressure (hPa)
   real(r8) pintmb(pcols,pverp)    ! Model interface pressure (hPa)

   real(r8) taua_lw(pcols,pver,nbndlw)     ! aerosol optical depth by band

   ! Set molecular weight ratios (for converting mmr to vmr)
   !  e.g. h2ovmr = h2ommr * amdw)
   real(r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
   real(r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
   real(r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
   real(r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
   real(r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
   real(r8), parameter :: amdo2 = 0.905140_r8   ! Molecular weight of dry air / oxygen
   real(r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
   real(r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

   real(r8), parameter :: dps = 1._r8/86400._r8 ! Inverse of seconds per day

   ! Cloud arrays for McICA 
   integer, parameter :: nsubclw = ngptlw       ! rrtmg_lw g-point (quadrature point) dimension
   integer :: permuteseed                       ! permute seed for sub-column generator

   real(r8) :: cicewp(pcols,pver)   ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,pver)   ! in-cloud cloud liquid water path
   real(r8) :: rei(pcols,pver)      ! ice particle effective radius (microns)
   real(r8) :: rel(pcols,pver)      ! liquid particle radius (micron)

   real(r8) cld_stolw(nsubclw, pcols, pver)     ! cloud fraction (mcica)
   real(r8) cicewp_stolw(nsubclw, pcols, pver)  ! cloud ice water path (mcica)
   real(r8) cliqwp_stolw(nsubclw, pcols, pver)  ! cloud liquid water path (mcica)
   real(r8) rei_stolw(pcols,pver)               ! ice particle size (mcica)
   real(r8) rel_stolw(pcols,pver)               ! liquid particle size (mcica)
   real(r8) tauc_stolw(nsubclw, pcols, pver)    ! cloud optical depth (mcica - optional)

   ! Includes extra layer above model top
   real(r8) uflx(pcols,pverp+1)  ! Total upwards longwave flux
   real(r8) uflxc(pcols,pverp+1) ! Clear sky upwards longwave flux
   real(r8) dflx(pcols,pverp+1)  ! Total downwards longwave flux
   real(r8) dflxc(pcols,pverp+1) ! Clear sky downwards longwv flux
   real(r8) hr(pcols,pver+1)     ! Longwave heating rate (K/d)
   real(r8) hrc(pcols,pver+1)    ! Clear sky longwave heating rate (K/d)
   !-----------------------------------------------------------------------

   ! mji/rrtmg

   ! Calculate cloud optical properties here if using CAM method, or if using one of the
   ! methods in RRTMG_LW, then pass in cloud physical properties and zero out cloud optical 
   ! properties here
   
   ! Zero optional cloud optical depth input array tauc_lw, 
   ! if inputting cloud physical properties into RRTMG_LW
   !          tauc_lw(:,:,:) = 0.
   ! Or, pass in CAM cloud longwave optical depth to RRTMG_LW
   ! do nbnd = 1, nbndlw
   !    tauc_lw(nbnd,:ncol,:pver) = cldtau(:ncol,:pver)
   ! end do

   ! Call mcica sub-column generator for RRTMG_LW

   ! Call sub-column generator for McICA in radiation
   call t_startf('mcica_subcol_lw')

   ! Select cloud overlap approach (1=random, 2=maximum-random, 3=maximum)
   icld = 2
   ! Set permute seed (must be offset between LW and SW by at least 140 to insure 
   ! effective randomization)
   permuteseed = 150

   ! These fields are no longer supplied by CAM.
   cicewp = 0.0_r8
   cliqwp = 0.0_r8
   rei = 0.0_r8
   rel = 0.0_r8

   call mcica_subcol_lw(lchnk, ncol, pver, icld, permuteseed, pmid, &
      cld, cicewp, cliqwp, rei, rel, tauc_lw, &
      cld_stolw, cicewp_stolw, cliqwp_stolw, rei_stolw, rel_stolw, tauc_stolw)

   call t_stopf('mcica_subcol_lw')

   
   call t_startf('rrtmg_lw')

   !
   ! Call RRTMG_LW model
   !
   ! Set input flags for cloud parameterizations
   ! Use separate specification of ice and liquid cloud optical depth.
   ! Use either Ebert and Curry ice parameterization (iceflglw = 0 or 1), 
   ! or use Key (Streamer) approach (iceflglw = 2), or use Fu method
   ! (iceflglw = 3), and Hu/Stamnes for liquid (liqflglw = 1).
   ! For use in Fu method (iceflglw = 3), rei is converted in RRTMG_LW
   ! from effective radius to generalized effective size using the
   ! conversion of D. Mitchell, JAS, 2002.  For ice particles outside
   ! the effective range of either the Key or Fu approaches, the 
   ! Ebert and Curry method is applied. 

   ! Input CAM cloud optical depth directly
   inflglw = 0
   iceflglw = 0
   liqflglw = 0
   ! Use E&C approach for ice to mimic CAM3
   !   inflglw = 2
   !   iceflglw = 1
   !   liqflglw = 1
   ! Use merged Fu and E&C params for ice
   !   inflglw = 2
   !   iceflglw = 3
   !   liqflglw = 1

   ! Convert incoming water amounts from specific humidity to vmr as needed;
   ! Convert other incoming molecular amounts from mmr to vmr as needed;
   ! Convert pressures from Pa to hPa;
   ! Set surface emissivity to 1.0 here, this is treated in land surface model;
   ! Set surface temperature
   ! Set aerosol optical depth to zero for now

   do k = 1, pver
      do i = 1, ncol
         h2ovmr(i,k) = (qnm(i,k) / (1._r8 - qnm(i,k))) * amdw
         o3vmr(i,k) = o3(i,k) * amdo
         co2vmr(i,k) = co2mmr(i,k) * amdc
         ch4vmr(i,k) = ch4(i,k) * amdm
         o2vmr(i,k)  = o2(i,k) * amdo2
         n2ovmr(i,k) = n2o(i,k) * amdn
         cfc11vmr(i,k) = cfc11(i,k) * amdc1
         cfc12vmr(i,k) = cfc12(i,k) * amdc2
         cfc22vmr(i,k) = 0._r8
         ccl4vmr(i,k) = 0._r8
         pmidmb(i,k) = pmid(i,k) * 1.e-2_r8
         pintmb(i,k) = pint(i,k) * 1.e-2_r8
         taua_lw(i,k,:nbndlw) = aer_lw_abs(i,k,:nbndlw)
      enddo
   enddo
   do i = 1, ncol
      emis(i,:nbndlw) = 1._r8
      tsfc(i) = tint(i,pverp)
      pintmb(i,pverp) = pint(i,pverp) * 1.e-2_r8
   enddo

   call rrtmg_lw(lchnk  ,ncol ,pver    ,icld    ,                 &
             pmidmb  ,pintmb  ,tnm     ,tint    ,tsfc    ,h2ovmr, &
             o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  ,cfc11vmr,cfc12vmr, &
             cfc22vmr,ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
             cld_stolw,tauc_stolw,cicewp_stolw,cliqwp_stolw ,rei, rel, &
             taua_lw, &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc)

   !
   !----------------------------------------------------------------------
   ! All longitudes: store history tape quantities
   ! Flux units are in W/m2 on output from rrtmg_lw and contain output for
   ! extra layer above model top with vertical indexing from bottom to top.
   ! Heating units are in K/d on output from RRTMG and contain output for
   ! extra layer above model top with vertical indexing from bottom to top.
   ! Heating units are converted to J/kg/s below for use in CAM. 
   !
   ! Reverse vertical indexing here for CAM arrays to go from top to bottom.
   !
   do i=1,ncol
      flwds(i) = dflx (i,1)
      fldsc(i) = dflxc (i,1)
      flns(i)  = uflx (i,1) - dflx (i,1)
      flnsc(i) = uflxc(i,1) - dflxc(i,1)
      flnt(i)  = uflx (i,pverp) - dflx (i,pverp)
      flntc(i) = uflxc(i,pverp) - dflxc(i,pverp)
      flut(i)  = uflx (i,pverp)
      flutc(i) = uflxc(i,pverp)
   end do
   do k = 1,pverp
      do i=1,ncol
         ful(i,k) = uflx(i,pverp+1-k)
         fdl(i,k) = dflx(i,pverp+1-k)
         fsul(i,k) = uflxc(i,pverp+1-k)
         fsdl(i,k) = dflxc(i,pverp+1-k)
      end do
   end do

   if (single_column.and.scm_crm_mode) then
      call outfld('FUL     ',ful,pcols,lchnk)
      call outfld('FDL     ',fdl,pcols,lchnk)
      call outfld('FULC    ',fsul,pcols,lchnk)
      call outfld('FDLC    ',fsdl,pcols,lchnk)
   endif

   do k = 1,pverp
      do i = 1,ncol
         fnl(i,k) = ful(i,k) - fdl(i,k)
         ! mji/ cam excluded this?
         fcnl(i,k) = fsul(i,k) - fsdl(i,k)
         !!
      end do
   end do

   ! Pass longwave heating to CAM arrays and convert from K/d to J/kg/s
   do k=1,pver
      do i=1,ncol
         qrl(i,k) = hr(i,pver+1-k)*cpair*dps
         qrlc(i,k) = hrc(i,pver+1-k)*cpair*dps
      end do
   end do
   ! Return 0 above solution domain
   if ( ntoplw > 1 )then
      qrl(:ncol,:ntoplw-1) = 0._r8
      qrlc(:ncol,:ntoplw-1) = 0._r8
   end if

   call t_stopf('rrtmg_lw')

end subroutine rad_rrtmg_lw

!-------------------------------------------------------------------------------

subroutine radlw_init()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize various constants for radiation scheme.
!
!-----------------------------------------------------------------------

   use hycoef, only : hypm

   integer :: k

   ! If the top model level is above ~90 km (0.1 Pa), set the top level to compute
   ! longwave cooling to about 80 km (1 Pa)
   if (hypm(1) .lt. 0.1_r8) then
      do k = 1, pver
         if (hypm(k) .lt. 1._r8) ntoplw  = k
      end do
   else
      ntoplw  = 1
   end if
   if (masterproc) then
      write(iulog,*) 'radlw_init: ntoplw =',ntoplw
   endif

   call rrtmg_lw_ini

end subroutine radlw_init

!-------------------------------------------------------------------------------

end module radlw
