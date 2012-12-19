c $Id$
c=======================================================================
!---! Energy-conserving sea ice model
!---! Code to compute ocean-ice basal flux and lateral melt and 
!---! These files bear a strong resemblance to those in CSIM4
!---! ice_therm_driver.F but unfortunately had to be modified a lot
!---! for this special use.
!---!
!---! author C. M. Bitz
!---!
!---! See Bitz, C.M., and W.H. Lipscomb, 1999: 
!---! An energy-conserving thermodynamic model of sea ice,
!---! J. Geophys. Res., 104, 15,669-15,677. 
!---!     
!---! The author grants permission to the public to copy and use this
!---! software without charge, provided that this Notice and any statement
!---! of authorship are reproduced on all copies and any publications that
!---! result from the use of this software must (1) refer to the publications 
!---! listed above and (2) acknowledge the origin and author of the model.
!---! This software is without warranty, expressed or implied, and the
!---! author assumes no liability or responsibility for its use. 
c=======================================================================
      module ice_ocn_flux

      use ice_kinds_mod
      use ice_constants, only: puny, c0, c1, c2, p5, Tfrez, saltz, rLfs,
     $                         rLfi, rcpi, rcpidepressT, rLfidepressT,
     $                         rhow, ni, cp_ocn, hi_min, Tffresh, pi,
     $                         salnew, plevmx
      use ice_dh,        only: energ
      
      implicit none

c=======================================================================

      contains

c=======================================================================

      subroutine init_frzmlt (sst, Tair, frzmlt, aice, hice, hs,
     &                        Tsfc, tiz, sal1d, Fbot, frazil, Rside, dt)

!---!-------------------------------------------------------------------
!---! initialize ocean-ice heat fluxes  bottom and lateral
!---! Assuming frzmlt is per grid box area
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::
     &   sst         ! SST (C)
     &,  Tair        ! Air temperature (K)
     &,  frzmlt      ! freeze melt potential from ocean (W/m**2)

      real (kind=dbl_kind), intent(inout) ::
     &   Tsfc        ! surface temperature (C)
     &,  Tiz(0:plevmx) ! local 1D ice temperature profile (C)
     &,  sal1d   (plevmx+1)     ! ice salinity (ppt)
     &,  aice        ! ice concentration
     &,  hice        ! ice thickness         (m)
     &,  hs          ! snow thickness (m) ! changed to intent(inout) OCT 11

      real (kind=dbl_kind), intent(out) ::
     &   Fbot        ! heat flx to ice bottom,         (W/m**2)
     &,  frazil      ! new ice growth rate
     &,  Rside       ! fraction of ice that melts from side
      real  (kind=dbl_kind), intent(in) :: dt   ! timestep

      real (kind=dbl_kind) ::
     &   Fsid        ! heat flx to side of ice (neg)  (W/m**2)
     &,  frcbot      ! fraction of available melt heat for bottom
     &,  frcsid      ! fraction of available melt heat for side
     &,  Fbotmx      ! maximum heat available for bottom melt (W/m**2)
     &,  Fsidmx      ! maximum heat available for side melt (W/m**2)
     &,  deltaT      ! amount above freezing temperature (C)
     &,  qice(ni)    ! ice energy of melting per unit volume (J/m**3)

      ! approx ustar from Steele 1995
      real (kind=dbl_kind), parameter ::
     &  ustar = 0.01_dbl_kind      ! skin friction velocity for Fbot (m/s)

      real (kind=dbl_kind), parameter ::
     &  Rfactor =  0.68_dbl_kind     ! shortwave penetration factor
     &,  zeta1  =  1.20_dbl_kind     ! shortwave e-folding depth
     &,  zeta2  = 28.00_dbl_kind     ! shortwave e-folding depth

      integer :: layer

      real (kind=dbl_kind) ::
     &    vii         ! initial volume of ice    (m)
     &,   vsn         ! initial volume of sno    (m)   ! ADDED OCT 11
     &,   vio         ! volume of new ice        (m)
     &,   vice        ! volume of old and new ice (m)
     &,   hnew        ! new ice thickness        (m)
     &,   Tnew        ! new ice Temperature (C)
     &,   ai0         ! new ice fraction
     &,   qio         ! new ice energy of melting per unit volume (J/m**3)
     &,   etot        ! total energy of melting (J/m**2)

      ! Parameters for basal and lateral heat flx 
      ! 0.006 = unitless param for basal heat flx ala McPhee and Maykut
      real (kind=dbl_kind), parameter :: 
     &   cpchr =-cp_ocn*rhow*0.006_dbl_kind

      ! params for lateral heat flx ala Maykut and Steel
      real (kind=dbl_kind), parameter :: 
     &   m1 = 1.6e-6_dbl_kind       ! (m/s/deg**m2)
     &,  m2 = 1.36_dbl_kind         ! unitless
     &,  floe_size = 300.0_dbl_kind       ! floe diameter (m)

      real (kind=dbl_kind) :: pi_eta

      pi_eta = pi/0.66_dbl_kind ! unitless 

      Rside=0._dbl_kind
!JR Init heat exchange between ice and underlying ocean
      Fbot = 0._dbl_kind
      frazil=0._dbl_kind

      if (frzmlt .ge. 0._dbl_kind) then     
      !-----------------------------------------------------------------
      ! freezing conditions, grow frazil ice 
      !-----------------------------------------------------------------
      ! so the ocean keep its temperature at or above freezing
      ! due to heat loss over open water
      !-----------------------------------------------------------------
           
        vii=aice*hice           ! initial ice volume
        vsn=aice*hs             ! initial sno volume   ! ADDED OCT 11

                                ! assume new ice has a energy of ice 
                                ! with 4psu at -1.8C
        qio = energ( Tfrez,salnew) 
        vio = -frzmlt*dt/qio    ! frazil ice volume

        frazil = vio/dt         ! frazil ice growth rate for history
        
        vice = vii+vio          ! combined ice volume
            
        ai0 = 1-aice            ! initial open water fraction
        if ( ai0 .gt. puny ) then
          hnew = vio/ai0        ! frazil ice thickness
          hnew = max(hi_min,hnew) ! make it >= hi_min   (do not skip this!)
          Tnew = min(sst,p5*(Tair-Tffresh)) ! Temp of frazil ice in C
          
          aice = aice + vio/hnew ! combine concentration of old ice and frazil
!cmb out sep4              Tsfc = Tsfc*aice +Tnew*ai0 ! combined ice surface temperature 
        else                    ! ai0 < puny
                                ! ocean wants to grow frazil ice but there is no open water
                                ! just increase the ice thickness 
          aice = 1.0_dbl_kind 
        endif
        hice = vice/aice        ! combined thickness

        hs = vsn/aice           ! ADDED OCT 11
                         
                                ! adjust the qi profile to reflect the increase in 
                                ! ice by frazil growth, add an equal amount to each layer
                                ! regardless of open water fraction
        if (vii.gt.puny) then
          do layer = 1,ni
            qice(layer) = energ(Tiz(layer),sal1d(layer))
            qice(layer) = (qice(layer)*vii+qio*vio)/vice
          end do
        else
          do layer = 1,ni
            qice(layer) = qio
          end do
        endif
!JR Changed qi (no declaration) to qice because that I am guessing that is what is needed
        call init_Tprofile(qice,tiz) ! recompute Tiz


      else

      !-----------------------------------------------------------------
      ! melting conditions
      !-----------------------------------------------------------------
      ! Compute fraction of heat available to melt side (frcsid)
      ! and that available to melt bottom (frcbot), and limit the
      ! amounts available according to frcsid*frzmlt and frcbot*frzmlt
      ! respectively (according to Briegleb I think)
      !-----------------------------------------------------------------

        if (aice.gt.c0) then

          frcbot = Rfactor*exp(-hice/zeta1) + 
     $            (c1-Rfactor)*exp(-hice/zeta2)
          frcsid = c1 - frcbot
          Fbotmx = frcbot*frzmlt
          Fsidmx = frcsid*frzmlt
          
            ! use boundary layer theory for Fbot
          deltaT = max((sst-Tfrez),c0) 
          Fbot = cpchr*deltaT*ustar ! < 0
          Fbot = max(Fbot,Fbotmx) ! Fbotmx < Fbot < 0

            ! allow for side melt
          etot = 0._dbl_kind
          do layer=1,ni
            etot = etot + energ(tiz(layer),saltz(layer))  
          enddo
!JR changed hi to hice because the variable hi does not exist.  This is a pure GUESS
!            etot = -hi*etot/ni + hs*rLfs  ! positive
          etot = -hice*etot/ni + hs*rLfs ! positive

          Fsid   = -etot*pi_eta/floe_size*m1*deltaT**m2
          Fsid   = max(Fsid,Fsidmx) ! Fsidmx < Fsid < 0
          Rside  = -Fsid*dt/etot
          Rside  = min(Rside,c1) ! positive

        endif                   ! aice
      endif                     ! sgn(frzmlt)
      
      end subroutine init_frzmlt

c=======================================================================

      subroutine lateral_melt (Rside, aice, hice, hs, qice,
     &                                Focn, meltl, dt)

!---!-------------------------------------------------------------------
!---! frazil ice growth / lateral melt
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::
     &   Rside       ! fraction of ice that melts from side

      real (kind=dbl_kind), intent(inout) ::
     &   aice        ! ice concentration
     &,  hice        ! ice thickness         (m)
     &,  hs          ! snow thickness (m)
!JR Added qice input arg declaration
     &,  qice(ni)    ! ice energy of melting per unit volume (J/m**3)
     &,  Focn        ! heat used for basal and lateral melt (W/m**2)

      real (kind=dbl_kind), intent(out) ::
     &   meltl       ! lateral melt rate diagnostic (m)

      real (kind=dbl_kind), intent(in) :: dt    ! model timestep

      integer :: layer
      real (kind=dbl_kind) ::
     &    vii         ! initial volume of ice    (m)
     &,   etot        ! total energy of melting (J/m**2)

      meltl = 0._dbl_kind

      if (aice.gt.puny .and.  Rside.gt.puny) then 

      !-----------------------------------------------------------------
      ! melt laterally (Rside >= 0)
      !-----------------------------------------------------------------
      ! etot=total energy of melt per unit area
        vii = hice*aice
        etot = 0._dbl_kind          
        do layer = 1,ni
          etot = etot + qice(layer)
        end do
        etot = -etot*vii/ni + aice*hs*rLfs ! positive

        Focn = Focn-Rside*etot/dt ! flux of heat used to melt from side 
        meltl = -Rside*vii/dt    ! lateral ice melt rate for history
        aice = aice*(c1-Rside)
        
      endif 
      
      end subroutine lateral_melt


c=======================================================================

      subroutine init_Tprofile(qi,tiz)

!---!-------------------------------------------------------------------
!---! compute the vertical temperature profile
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) :: 
     &   qi(ni)  ! enthalpy per unit volume     (J/m**3)

      real (kind=dbl_kind), intent(out) :: 
     &  tiz(0:ni)  ! temp of each layer    (C)

      integer :: layer

      real (kind=dbl_kind) ::
     &      T2        ! solution of the quadratic eq that gets used (C)
     &,     q, B, C   ! variables to help solve quad. eq.
     &,     B_2, root

      ! compute the midpoint temperature from the 
      ! energy of melting for each layer
      ! by solving the quadratic equation

!JR Changed loop limit to ni because nilay(nc) is not available
      do layer = 1,ni
        q = qi(layer)+rLfi-rcpidepressT*saltz(layer)
        B = -q/rcpi
        C = -rLfidepressT*saltz(layer)/rcpi
        B_2  = B/c2
        root = sqrt(B_2*B_2-C)
c        T1   = -B_2+root
        T2   = -B_2 - root
        tiz(layer) = T2
      enddo

      end subroutine init_Tprofile

c=======================================================================

      end module ice_ocn_flux

c=======================================================================
