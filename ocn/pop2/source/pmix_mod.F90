!  CVS: $Id: pmix_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module pmix_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!     Calculating Richardson number riu (at U/V-point);
!                                   rit (at T-point)
!     in order to computing vertical mixing coefficients based on
!     Pacanowski & Philander (JPO vol 11, #11, 1981).
!
!     Note: this parameterization was designed for equatorial models
!     and may not do a good job in mid or high latitudes. Simulations
!     in these regions (where vertical shear is small) are improved with
!     the addition of solar short wave penetration into the ocean which
!     reduces buoyancy and enhances vertical mixing.
!
      INTEGER:: RTST,RTEND,RUST,RUEND
!
      real(r8):: wndmix, fricmx, diff_cbt_back, diff_cbt_limit,&
                     visc_cbu_back, visc_cbu_limit
!lhl1204 real(r8),dimension(:,:,:),allocatable:: ric,rit
      real(r8),dimension(:,:,:),allocatable:: ric,rict
      real(r8),dimension(:,:,:),allocatable:: rit,riu
!      real(r8),dimension(imt,jmt,kmm1):: ricdt,ricdu
      real(r8),dimension(imt,jmt,kmm1):: ricdt
!      real(r8),dimension(imt,jmt,kmm1):: ridu,ridt,s2u,s2t
      real(r8),dimension(imt,jmt,kmm1):: ridt,s2u,s2t
!lhl1204
!
!
!
!
!
!
!-----------------------------------------------------------------------
!     Solar Shortwave energy penetrates below the ocean surface. Clear
!     water assumes energy partitions between two exponentials as
!     follows:
!
!     58% of the energy decays with a 35 cm e-folding scale
!     42% of the energy decays with a 23 m e-folding scale
!
!     if the thickness of the first ocean level "dzt(1)" is 50 meters,
!     then shortwave penetration wouldn't matter. however, for
!     dzt(1) = 10 meters, the effect can be significant. this can be
!     particularly noticable in the summer hemisphere.
!
!     Paulson and Simpson (1977 Irradiance measurements in the upper
!                               ocean JPO 7, 952-956)
!     Also see ... Jerlov (1968 Optical oceanography. Elsevier)
!                  A General Circulation Model for Upper Ocean
!                  Simulaton (Rosati and Miyakoda JPO vol 18,Nov 1988)
!   based on chlorophyll concentration shortwave radiation
!   can see Ohlmann Ocean Radiant Heating in Climate models( JC 5 2003 )
!-----------------------------------------------------------------------
!
#if (defined SOLAR)
      REAL(r8):: pen(kmm1)
#endif

#if (defined SOLARCHLORO)
      real(r8),dimension(imt,jmt,km):: pen_chl   !  be different from pen(kmm1)
#endif

end module pmix_mod
