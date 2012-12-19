!  CVS: $Id: yy00.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ===============
      SUBROUTINE YY00
!     ===============
 
#include <def-undef.h>
use param_mod
use output_mod

      IMPLICIT NONE

!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = 1,JMT
         DO I = 1,IMT
            Z0MON (I,J)= 0.0
            HIMON (I,J)= 0.0
            HDMON (I,J)= 0.0
            ICMON (I,J,1)= 0.0
            ICMON (I,J,2)= 0.0
         END DO
 
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               TSMON (I,J,K)= 0.0
               SSMON (I,J,K)= 0.0
               USMON (I,J,K)= 0.0
               VSMON (I,J,K)= 0.0
               WSMON (I,J,K)= 0.0
#if (defined SMAG_OUT)
               AM3MON (I,J,K)= 0.0
#endif
            END DO
         END DO
      END DO

!$OMP PARALLEL DO PRIVATE (N,K,J,I)
      DO N = 1,NTRA
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               trendmon (I,J,K,N)= 0.0
               axmon (I,J,K,N)= 0.0
               aymon (I,J,K,N)= 0.0
               azmon (I,J,K,N)= 0.0
               dxmon (I,J,K,N)= 0.0
               dymon (I,J,K,N)= 0.0
               dzmon (I,J,K,N)= 0.0
!
               ddymon (I,J,K,N)= 0.0
#ifdef ISO
               axmon_iso (I,J,K,N)= 0.0
               aymon_iso (I,J,K,N)= 0.0
               azmon_iso (I,J,K,N)= 0.0
               dxmon_iso (I,J,K,N)= 0.0
               dymon_iso (I,J,K,N)= 0.0
               dzmon_iso (I,J,K,N)= 0.0
!
               aaymon_iso (I,J,K,N)= 0.0
               ddymon_iso (I,J,K,N)= 0.0
#endif
            END DO
         END DO
      END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (N,J,I)
      DO N = 1,NTRA
         DO J = 1,JMT
            DO I = 1,IMT
               netmon (I,J,N)= 0.0
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               penmon (I,J,K)= 0.0
            END DO
         END DO
      END DO

!lhl1204
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = 1,JMT
            DO I = 1,IMT
               mldmon (I,J)= 0.0
            END DO
         END DO
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               akmmon (I,J,K)= 0.0
               aktmon (I,J,K)= 0.0
               aksmon (I,J,K)= 0.0
            END DO
         END DO
      END DO

!lhl1204
 
      RETURN
      END SUBROUTINE YY00
 
 
