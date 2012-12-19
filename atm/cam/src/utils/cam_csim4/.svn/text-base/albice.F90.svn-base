
subroutine albice(lchnk   ,ncol    ,Tair    ,snowh          ,&
                  asdir   ,aldir   ,asdif   ,aldif, icefrac ,&
                  sicthk )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute surface albedos
!
! Method: 
! Computes surface albedos for direct/diffuse incident radiation for
! two spectral intervals:
!   s = 0.2-0.7 micro-meters
!   l = 0.7-5.0 micro-meters
!
! Albedos specified as follows:
! Ocean with      Surface albs specified; combined with overlying snow
!   sea ice       
!
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols
  use ice_constants,only: c0, c1, p2, p33, rhofresh, rhos, Tffresh, snowpatch, timelt
  use dycore, only: dycore_is, get_resolution

  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: lchnk            ! chunk identifier
  integer , intent(in) :: ncol             ! number of atmospheric columns

  real(r8), intent(in) :: Tair(pcols)      ! bottom level air temp
  real(r8), intent(in) :: snowh(pcols)     ! Snow depth (liquid water equivalent)
  real(r8), intent(in) :: icefrac(pcols)   ! Areal sea-ice fraction
  real(r8), intent(in) :: sicthk(pcols)    ! sea-ice thickness
  real(r8), intent(out):: asdir(pcols)     ! Srf alb for direct rad   0.2-0.7 micro-ms
  real(r8), intent(out):: aldir(pcols)     ! Srf alb for direct rad   0.7-5.0 micro-ms
  real(r8), intent(out):: asdif(pcols)     ! Srf alb for diffuse rad  0.2-0.7 micro-ms
  real(r8), intent(out):: aldif(pcols)     ! Srf alb for diffuse rad  0.7-5.0 micro-ms

!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                  ! Longitude index
                             ! albedos for ice in each category
  real(r8) :: alvdrn (pcols) ! visible, direct   (fraction)
  real(r8) :: alidrn (pcols) ! near-ir, direct   (fraction)
  real(r8) :: alvdfn (pcols) ! visible, diffuse  (fraction)
  real(r8) :: alidfn (pcols) ! near-ir, diffuse  (fraction)
!-----------------------------------------------------------------------
  real (r8), parameter :: albocn = 0.06_r8  ! ocean albedo
  real (r8), parameter :: &
       ahmax    = 1.0_r8,   &! thickns above which ice alb is const,
       albicev  = 0.68_r8,  &! visible ice albedo for h > ahmax
       albicei  = 0.30_r8,  &! near-ir ice albedo for h > ahmax
       dT_mlt   = 1._r8,    &! change in temp to give dalb_mlt change
       dalb_mlt = -0.075_r8,&! albedo change per dT_mlt change
       dalb_mltv= -0.100_r8,&! albedo vis change per dT_mlt change in temp for snow
       dalb_mlti= -0.150_r8  ! albedo nir change per dT_mlt change in temp for snow

  real (r8) :: &
       albsnowv,            &! cold snow albedo, visible
       albsnowi              ! cold snow albedo, near IR

! parameter for fractional snow area 
  real(r8)  fhtan ! factor used in dependence of albedo on ice thickness
  real(r8)  vicen(pcols),vsnon(pcols),aicen(pcols),tsfcn(pcols)
  real (r8) hi    ! ice thickness  (m)
  real (r8) hs    ! snow thickness (m)
  real (r8) snw   !
  real (r8) albo  ! effective ocean albedo, function of ice thickness
  real (r8) asnow ! snow-covered area fraction
  real (r8) asnwv ! snow albedo, visible 
  real (r8) asnwi ! snow albedo, near IR
  real (r8) fh    ! piecewise linear function of thickness 
  real (r8) fT    ! piecewise linear function of surface temperature
  real (r8) dTs   ! difference of Tsfc and Timelt
!-----------------------------------------------------------------------

!
! Set snow albedos (every time, since there is no initialization procedure)
!
  if ( dycore_is ('LR') ) then
     albsnowv = 0.96_r8 ! cold snow albedo, visible
     albsnowi = 0.68_r8 ! cold snow albedo, near IR
  else 
     if ( get_resolution() == 'T31' ) then
        albsnowv = 0.91_r8 ! cold snow albedo, visible
        albsnowi = 0.63_r8 ! cold snow albedo, near IR
     else
        albsnowv = 0.96_r8 ! cold snow albedo, visible
        albsnowi = 0.68_r8 ! cold snow albedo, near IR
     endif
  endif
  
!
! Initialize all sea ice surface albedos to zero
!
  asdir(:) = 0._r8
  aldir(:) = 0._r8
  asdif(:) = 0._r8
  aldif(:) = 0._r8
  alvdrn(:) = 0._r8
  alidrn(:) = 0._r8
  alvdfn(:) = 0._r8
  alidfn(:) = 0._r8

  fhtan = atan(ahmax*5._r8) 

  do i=1,ncol
     if (icefrac(i) > 0._r8) then
        hi  = sicthk(i)
        snw = snowh(i)*rhofresh/rhos
        aicen(i) = icefrac(i)
        vicen(i) = hi*aicen(i) 
        !---------------------------------------------------------
        ! keep snow/ice boundary above sea level by reducing snow
        !---------------------------------------------------------
        vsnon(i) = min(snw*aicen(i),p33*vicen(i))
        Tsfcn(i) = min(Tair(i)-Tffresh,-p2)   ! deg C       
        !---------------------------------------------------------
        ! make linear temp profile and compute enthalpy
        !---------------------------------------------------------

        hi = vicen(i) / aicen(i)
        hs = vsnon(i) / aicen(i)
        
        ! bare ice, thickness dependence
        fh = min(atan(hi*5._r8)/fhtan,c1)
        albo = albocn*(c1-fh)
        alvdfn(i) = albicev*fh + albo
        alidfn(i) = albicei*fh + albo

        ! bare ice, temperature dependence
        dTs = Timelt - Tsfcn(i)
        fT = min(dTs/dT_mlt-c1,c0)
        alvdfn(i) = alvdfn(i) - dalb_mlt*fT
        alidfn(i) = alidfn(i) - dalb_mlt*fT

        if( hs .gt. 0._r8 ) then
           ! fractional area of snow on ice (thickness dependent)
           asnow = hs / ( hs + snowpatch ) 
           asnwv = albsnowv
           asnwi = albsnowi
           ! snow on ice, temperature dependence
           asnwv = asnwv - dalb_mltv*fT
           asnwi = asnwi - dalb_mlti*fT
           
           ! combine ice and snow albedos
           alvdfn(i) = alvdfn(i)*(c1-asnow) + asnwv*asnow
           alidfn(i) = alidfn(i)*(c1-asnow) + asnwi*asnow
        endif
        alvdrn(i) = alvdfn(i)
        alidrn(i) = alidfn(i)
     endif  ! aicen > 0._r8

     if (icefrac(i) > 0._r8) then
        asdir(i)  = alvdrn(i)
        aldir(i)  = alidrn(i)
        asdif(i) = alvdfn(i)
        aldif(i) = alidfn(i)
     end if
  enddo
!
  return
end subroutine albice

