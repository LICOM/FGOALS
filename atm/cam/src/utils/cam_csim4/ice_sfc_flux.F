c=======================================================================
!---! atmospheric boundary interface (stability based flux calculations)
!---!
!---! based on code by Bill Large
c=======================================================================

      module ice_sfc_flux

      use ice_kinds_mod
      use ice_constants, only: cp_air, cpwv, c0, c1, c2, c4, c5, c8, 
     &                         c10, c16, Tffresh, iceruf, zref, vonkar,
     &                          Lsub, stefan_boltzmann, gravit, p5
      use shr_const_mod, only: SHR_CONST_ZVIR, SHR_CONST_RDAIR
      use ppgrid,        only: pcols

      implicit none

c=======================================================================

      contains

c=======================================================================

      subroutine ice_atm_flux (Tsf, uatm, vatm, tair,
     &                         qa    ,potT   ,zlvl    ,pbot  , Flw ,
     &                         swvdr   ,swidr   ,swvdf   ,swidf , 
     &                         alvdr   ,alidr   ,alvdf   ,alidf , 
     &                         qqq, TTT, emissivity,
     &                         strx, stry, flwup, fsh, flh, Tref,
     &                         flwdabs, fswabs, fswabsv, fswabsi,
     &                         dflhdT, dfshdT, dflwdT,
     &                         npts, indx)

!---!-------------------------------------------------------------------
!---! Compute ice (or ocn)-atm surface fluxes, stress, and reference 
!---! temperature
!---! NOTE: 
!---!  o all fluxes are positive downward
!---!  o net heat flux = fswabs + flwup + (1-emissivity)flwdn + fsh(i) + flh
!---!  o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!---!  o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
!---!
!---! ASSUME:
!---!  o The saturation humidity of air at T(K): qsat(T)  (kg/m**3)
!---!
!---! code originally based on CSM1
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in), dimension(pcols) ::
     &   Tsf    ! surface temperature of ice 
     &,  uatm  ! surface u wind
     &,  vatm  ! surface v wind
     &,  tair        ! Bottom level temperature
     &,  qa          ! Bottom level specific humidity
     &,  potT         ! Bottom level potential temperature
     &,  zlvl          ! Bottom level height above surface
     &,  pbot        ! Bottom level pressure
     &,  Flw     ! net down longwave radiation at surface
     &,  swvdr   ! direct beam solar radiation onto srf (sw)
     &,  swidr   ! direct beam solar radiation onto srf (lw)
     &,  swvdf   ! diffuse solar radiation onto srf (sw)
     &,  swidf   ! diffuse solar radiation onto srf (lw)
     &,  alvdr   ! ocean + ice albedo: shortwave, direct
     &,  alvdf   ! ocean + ice albedo: shortwave, diffuse
     &,  alidr   ! ocean + ice albedo: longwave, direct
     &,  alidf   ! ocean + ice albedo: longwave, diffuse
      real (kind=dbl_kind), intent(in) ::
     &   qqq     ! for qsat, dqsatdt
     &,  TTT     ! for qsat, dqsatdt
      real (kind=dbl_kind), intent(in) :: emissivity
      

      real (kind=dbl_kind), intent(out), dimension(pcols) ::
     &   strx   ! x surface stress (N)
     &,  stry   ! y surface stress (N)
     &,  Tref    ! reference height temperature  (K)
     &,    Flwdabs   ! down long-wave  absorbed heat flx   (W/m**2)
     &,    Flwup     ! emitted long-wave upward heat flux  (W/m**2)
     &,    fswabs     ! fswabs sum
     &,    fswabsv    ! fswabs(i) in vis (wvlngth < 700nm)  (W/m**2)
     &,    fswabsi    ! fswabs(i) in nir (wvlngth > 700nm)  (W/m**2)
     &,    fsh      ! sensible         heat flux  (W/m**2)
     &,    flh      ! latent           heat flux  (W/m**2)
     &,    dflhdT     ! d(flh(i))/d(T)      (W/m**2/K)
     &,    dfshdT     ! d(fsh(i))/d(T)      (W/m**2/K)
     &,    dFlwdT    ! d(Flwup(i))/d(T)     (W/m**2/K)
      integer, intent(in) :: npts
      integer, intent(in) :: indx(*)


! local variables
      real (kind=dbl_kind) :: 
     &   dssqdt ! derivative of ssq wrt Ti (kg/kg/K)
     &,  delt   ! potential T difference   (K)
     &,  delq   ! humidity difference      (kg/kg)
 
      integer :: k      ! iteration index
      integer :: i
      integer :: ii
      real (kind=dbl_kind) :: 
     &   TsfK   ! surface temperature in Kelvin (K)
     &,  thva   ! virtual temperature      (K)
     &,  stable ! stability factor
     &,  rdn    ! sqrt of neutral exchange coefficient (momentum)
     &,  rhn    ! sqrt of neutral exchange coefficient (heat)
     &,  ren    ! sqrt of neutral exchange coefficient (water)
     &,  hol    ! H (at zlvl(i)  ) over L
     &,  xsq    ! temporary variable
     &,  xqq    ! temporary variable
     &,  psimh  ! stability function at zlvl(i)   (momentum)
     &,  psixh  ! stability function at zlvl(i)   (heat and water)
     &,  alz    ! ln(zlvl(i)  /z10)
     &,  tau    ! stress at zlvl(i)
     &,  bn     ! exchange coef funct for interpolation
     &,  bh     ! exchange coef funct for interpolation
     &,  fac    ! interpolation factor
     &,  ln0    ! log factor for interpolation
     &,  ln3    ! log factor for interpolation
     &,  ustar  ! ustar (m/s)
     &,  tstar  ! tstar
     &,  qstar  ! qstar
     &,  rd     ! sqrt of exchange coefficient (momentum)
     &,  re     ! sqrt of exchange coefficient (water)            
     &,  rh     ! sqrt of exchange coefficient (heat)
     &,  vmag   ! surface wind magnitude   (m/s)
     &,  ssq    ! sat surface humidity     (kg/kg)
     &,  cp     ! specific heat of moist air
     &,  rhoa   ! air density (kg/m**3)
     &,  lhcoef
     &,  shcoef
     &,  wind  ! surface wind speed 

      real (kind=dbl_kind), parameter ::
     &   cpvir = cpwv/cp_air - c1  ! Defined as cpwv/cp_air - 1.
     &,  zTref  = c2          ! reference height for air temperature (m)
     &,  umin  = c1          ! minimum wind speed (m/s)
     &,  zvir  = SHR_CONST_ZVIR    ! rh2o/rair - 1.0

      ! local functions
      real (kind=dbl_kind) :: 
     &   Tk      ! temperature (K)
     &,  qsat    ! the saturation humididty of air (kg/m**3)
     &,  dqsatdt ! derivative of qsat wrt surface temperature
     &,  xd      ! dummy argument  
     &,  psimhu  ! unstable part of psimh
     &,  psixhu  ! unstable part of psimx

      qsat(Tk)    = qqq / exp(TTT/Tk)

      dqsatdt(Tk) = (TTT / Tk**2) * qqq / exp(TTT/Tk)

      psimhu(xd)  = log((c1+xd*(c2+xd))*(c1+xd*xd)/c8)
     $              - c2*atan(xd) + 1.571_dbl_kind

      psixhu(xd)  =  c2 * log((c1 + xd*xd)/c2)

      ! define some needed variables
!cdir nodep
      do ii = 1, npts
        i = indx(ii)
      TsfK   = Tsf(i) +Tffresh                 !  surface temp (K)
!JR      if (TsfK < 200. .or. TsfK > 350.) then
!JR        write(6,*)'ice_sfc_flux: bad TsfK=',TsfK
!JR        call endrun
!JR      end if
      wind   = sqrt(uatm(i)*uatm(i) + vatm(i)*vatm(i))
      vmag   = max(umin, wind)
      rhoa   = pbot(i)/(SHR_CONST_RDAIR*tair(i))
      thva   = potT(i) * (c1 + zvir * qa(i))      ! virtual pot temp (K)
      ssq    = qsat   (TsfK) / rhoa         ! sat surf hum (kg/kg)
      dssqdt = dqsatdt(TsfK) / rhoa         ! deriv of ssq wrt Ti 
      delt   = potT(i) - TsfK                  ! pot temp diff (K)
      delq   = qa(i) - ssq                     ! spec hum dif (kg/kg)
      alz    = log(zlvl(i)/zref) 
      cp     = cp_air*(c1 + cpvir*ssq)
ccc        write(6,*) 'cp is the problem',cp,ssq,qsat(TsfK),TsfK,rhoa, 
ccc       & Tsf(i),Tffresh
ccc  
ccc        write(6,*) 'IN ice_sfc_flux', Tsf(i), uatm(i), vatm(i), tair(i),
ccc       &    qa(i)    ,potT(i)   ,zlvl(i)    ,pbot(i)  , Flw(i) ,
ccc       &    swvdr(i)   ,swidr(i)   ,swvdf(i)   ,swidf(i) , 
ccc       &    alvdr(i)   ,alidr(i)   ,alvdf(i)   ,alidf(i) 
      !------------------------------------------------------------
      ! first estimate of Z/L and ustar, tstar and qstar
      !------------------------------------------------------------

      ! neutral coefficients, z/L = 0.0 
      rdn = vonkar/log(zref/iceruf)
      rhn = rdn
      ren = rdn

      ! ustar,tstar,qstar
      ustar = rdn * vmag
      tstar = rhn * delt  
      qstar = ren * delq  

      !------------------------------------------------------------
      ! iterate to converge on Z/L, ustar, tstar and qstar
      !------------------------------------------------------------

!cdir expand=5
      do k=1,5

        ! compute stability & evaluate all stability functions 
        hol    = vonkar * gravit * zlvl(i)
     $           * (tstar/thva+qstar/(c1/zvir+qa(i))) / ustar**2
        hol    = sign( min(abs(hol),c10), hol )
        stable = p5 + sign(p5 , hol)
        xsq    = max(sqrt(abs(c1 - c16*hol)) , c1)
        xqq    = sqrt(xsq)
        psimh  = -c5*hol*stable + (c1-stable)*psimhu(xqq)
        psixh  = -c5*hol*stable + (c1-stable)*psixhu(xqq)

        ! shift all coeffs to measurement height and stability
        rd = rdn / (c1+rdn/vonkar*(alz-psimh))
        rh = rhn / (c1+rhn/vonkar*(alz-psixh))
        re = ren / (c1+ren/vonkar*(alz-psixh))

        ! update ustar, tstar, qstar using updated, shifted coeffs 
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 

      enddo    ! end iteration

      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------

      shcoef = rhoa*ustar*cp  *rh
      lhcoef = rhoa*ustar*Lsub*re

      !------------------------------------------------------------
      ! momentum flux
      !------------------------------------------------------------
      ! tau = rhoa * ustar * ustar 
      ! strx(i) = tau * uatm(i) / vmag 
      ! stry(i) = tau * vatm(i) / vmag 
      !------------------------------------------------------------

      tau = rhoa * ustar * rd    ! not the stress at zlvl(i)
      strx(i) = tau * uatm(i) 
      stry(i) = tau * vatm(i)

      !------------------------------------------------------------
      ! reference temperature interpolation
      !------------------------------------------------------------
      ! Assume that 
      ! cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
      ! 1/sqrt(cn)=1/rdn and sqrt(cm)/ch=1/rh 
      !------------------------------------------------------------
      bn = vonkar/rdn
      bh = vonkar/rh

      ! Interpolation factor for stable and unstable cases
      ln0 = log(c1 + (zTref/zlvl(i))*(exp(bn) - c1))
      ln3 = log(c1 + (zTref/zlvl(i))*(exp(bn - bh) - c1))
      fac = (ln0 - zTref/zlvl(i)*(bn - bh))/bh 
     &         * stable
     $    + (ln0 - ln3)/bh * (c1-stable)
      fac = min(max(fac,c0),c1)

      Tref(i) = TsfK + (tair(i) - TsfK)*fac

      ! shortwave radiative flux
      fswabsv(i)  = swvdr(i)*(c1-alvdr(i)) + swvdf(i)*(c1-alvdf(i))
      fswabsi(i)  = swidr(i)*(c1-alidr(i)) + swidf(i)*(c1-alidf(i))
      fswabs(i)   = fswabsv(i) + fswabsi(i)

      ! longwave radiative flux
      Flwdabs(i) = emissivity*Flw(i)
      Flwup(i)   = -emissivity*stefan_boltzmann * TsfK**4
!JR      write(6,*)'ICE_ATM_FLUX: Flwup(i),TsfK=',Flwup(i),TsfK

      ! downward latent and sensible heat fluxes
      flh(i) = lhcoef * delq
      fsh(i) = shcoef * delt

!      write(6,*) '(ice_sfc_flux)',fsh(i),shcoef,delt

      ! derivatives wrt surface temp
      dFlwdT(i) = - emissivity*stefan_boltzmann * c4*TsfK**3 
      dflhdT(i) = - lhcoef * dssqdt
      dfshdT(i) = - shcoef

!! NOTE
! if the temperature dependence of flh(i), fsh(i) and Flwup(i) are
! included in the iteration for the ice sfc temperature, move
! fsw* out of this routine and use the following expressions instead:
! (NONE of the intent(in) variables will then be needed)

c**        Qsfc = Qcoef * exp(22.47_dbl_kind*(c1-Tffresh/TsfK))
c**        fsh(i)  = shcoef*(potT(i) - TsfK)
c**        flh(i)  = lhcoef*(qa(i) - Qsfc)

c**        dfshdT(i) = - shcoef
c**        dflhdT(i) = - lhcoef*lvrrv*Qsfc/TsfK**2
      end do ! ii = 1, npts

      end subroutine ice_atm_flux

c=======================================================================

      end module ice_sfc_flux

c=======================================================================

