!  CVS: $Id: invtri.F90,v 1.5 2003/08/12 09:06:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     =================
      SUBROUTINE INVTRI (WK,TOPBC,DCB,AIDIF,C2DTTS)
!     =================
use param_mod
use pconst_mod
      IMPLICIT NONE
      REAL    :: WK (IMT,JMT,KM),TOPBC (IMT,JMT),DCB (IMT,JMT,KM)
      REAL    :: A8 (KM),B8 (KM),C8 (KM),D8 (KM),E8 (0:KM),F8 (0:KM)
      REAL    :: AIDIF,C2DTTS,G0
      INTEGER :: KZ
 
!!$OMP PARALLEL DO PRIVATE (J,I,K,KZ,G0)
      JJJ: DO J = 2,JMM
         III: DO I = 2,IMM
            IF (ITNU (I,J) == 0) cycle III
            DO K = 2,KM  !lyc
               A8 (K) = DCB (I,J,K -1)* ODZT (K )* ODZP (K)* C2DTTS * AIDIF
               D8 (K) = WK (I,J,K)
            ENDDO
            DO K=2,KM-1
               C8 (K) = DCB (I,J,K )* ODZT (K+1 )* ODZP (K)* C2DTTS * AIDIF
               B8 (K) = 1.0+ A8 (K) + C8 (K)
               E8 (K -1) = 0.0
               F8 (K -1) = 0.0
            END DO
 
!     B. C. AT TOP
            K = 1
            A8 (K) = ODZP (K)* C2DTTS * AIDIF
            C8 (K) = DCB (I,J,K)* ODZT (K+1)* ODZP (K)* C2DTTS * AIDIF
            B8 (K) = 1.0+ C8 (K)
            D8 (K) = WK (I,J,K)
            E8 (K -1) = 0.0
            F8 (K -1) = 0.0
!     B. C. AT BOTTOM
            KZ = ITNU (I,J)
            IF (KZ /= 0) THEN
               B8 (KZ) = 1.0+ A8 (KZ)
               C8 (KZ) = ODZP (KZ)* C2DTTS * AIDIF
               E8 (KZ) = 0.0
               F8 (KZ) = 0.0
            END IF
 
!     NOW INVERT
            DO K = KM,1, -1
               IF (K <= ITNU (I,J)) THEN
                  G0 = 1.0/ (B8 (K) - C8 (K)* E8 (K))
                  E8 (K -1) = A8 (K)* G0
                  F8 (K -1) = (D8 (K) + C8 (K)* F8 (K))* G0
               END IF
 
            END DO
 
!     B.C. AT SURFACE
            WK (I,J,1) = (E8 (0)* TOPBC (I,J) + F8 (0))* VIT (I,J,1)
            DO K = 2,KM
               WK (I,J,K)= (E8 (K -1)* WK (I,J,K -1) + F8 (K -1))* VIT (I,J,K)
            END DO
 
         END DO III
      END DO JJJ
      RETURN
      END SUBROUTINE INVTRI
 
 
#else
      SUBROUTINE INVTRI ()
      RETURN
      END SUBROUTINE INVTRI
#endif 
