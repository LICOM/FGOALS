!  CVS: $Id: vinteg.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ==========================
      SUBROUTINE VINTEG (WK3,WK2)
!     ==========================
use precision_mod 
use param_mod
use pconst_mod
 
      IMPLICIT NONE
      REAL(r8):: WK3 (IMT,JMT,KM),WK2 (IMT,JMT)
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            WK2 (I,J)= 0.0D0
         END DO
      END DO
 
 
      DO K = 1,KM
!$OMP PARALLEL DO PRIVATE (J,I)
        DO J = JST,JET
         DO I = 1,IMT
               WK2 (I,J)= WK2 (I,J) + DZP (K)* OHBU(I,J)* WK3 (I,J,K)  &
                                  * VIV (I,J,K)
            END DO
         END DO
      END DO
 
 
      RETURN
      END SUBROUTINE VINTEG
 
 
