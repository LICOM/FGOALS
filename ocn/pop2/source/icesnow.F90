!  CVS: $Id: icesnow.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ==================
      SUBROUTINE ICESNOW
!     ==================
!     Sea Ice Model
#include <def-undef.h> 
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
      IMPLICIT NONE
!
      real(r8) :: sal_ocn, sal_ice, tdiff, heat_ice_fusion
      real(r8) ::  t_mix, s_mix
!
      heat_ice_fusion = 3.337d+5        !J/Kg
      sal_ocn=35.0D0                     !Reference salinity for ocean
      sal_ice= 0.0D0                     !Reference salinity for sea ice
!
!------------------------------------------------------------
!  if SST exceed -1.8C to restore it to -1.8C
!------------------------------------------------------------
 
 
      KKK : DO K = 1, 1
      JJJ : DO J = JST,JMT
         III : DO I = 1,IMT
 
            IF (ITNU (I,J) == 0) CYCLE III
 
            IF (AT (I,J,K,1) < TBICE) THEN
#ifdef  COUP
               tdiff        = TBICE- AT(i,j,k,1)
               licomqice (I,J  ) = licomqice(i,j)+tdiff*dzp(k)/dzp(1)
               at(i,j,1,2)=at(i,j,1,2)+tdiff*(sal_ocn-sal_ice)  &
                          *CP/heat_ice_fusion*0.001D0
#endif
               AT (I,J,k,1) = TBICE
            END IF
 
         END DO  III
      END DO JJJ
      END DO KKK
!
#ifdef  COUP
      DO J=JST,JMT
      DO I=1,IMT
         IF (licomqice(I,J) > 0.0 .and. at(i,j,1,1) > TBICE) THEN
              tdiff=min((at(i,j,1,1)-tbice),licomqice(i,j))
              licomqice(i,j)=licomqice(i,j)-tdiff
              at(i,j,1,1)=at(i,j,1,1)-tdiff
              at(i,j,1,2)=at(i,j,1,2)-tdiff*(sal_ocn-sal_ice)   &
                         *CP/heat_ice_fusion*0.001D0
          END IF
      END DO
      END DO 
#endif
!
      DO K=2,KM
      DO J=JST,JMT
      DO I=1,IMT
         IF (AT(I,J,K,1) < TBICE) AT(I,J,K,1)=TBICE
      ENDDO
      ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE ICESNOW
 
 
