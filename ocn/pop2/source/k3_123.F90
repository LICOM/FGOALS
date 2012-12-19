!  CVS: $Id: k3_123.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     =================
      SUBROUTINE K3_123
!     =================
 
!     compute K2(,,,1:3) at the center of the bottom face of "T" cells
!     use "c1e10" to keep the exponents in range.
use precision_mod 
use param_mod
use pconst_mod
use isopyc_mod
 
      IMPLICIT NONE
      REAL(r8) :: c1e10,eps,p5,p25,c0,c1,chkslp,ahfctr,xx
!lhl060506
      REAL(r8) , dimension(IMT,JMT,KM) :: F1,F2
      REAL(r8) :: NONDIMR,SLOPEMOD

      F1=0.0D0
      F2=0.0D0

!lhl060506
 
!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------
 
      c1e10 = 1.0d10
 
      eps = 1.0d-25
      p5 = 0.5D0
      p25 = 0.25D0
      c0 = 0.0D0
      c1 = 1.0D0
!yyq 0803408
!M
 
!$OMP PARALLEL DO PRIVATE (j,k,m,i)
      DO j = 2,jmm
         DO k = 2,km
            m = kisrpl (k)
            DO i = 2,imm
               e (i,k -1,j,1) = OTX (j)* p25* c1e10 &
               * (tmask (i -1,k -1,j)* tmask (i,k -1,j)* (rhoi (i,k -1,j,m) &
               - rhoi (i -1,k -1,j,m)) &
               + tmask (i,k -1,j)* tmask (i +1,k -1,j)* (rhoi (i +1,k -1,j,m)&
               - rhoi (i,k -1,j,m)) &
               + tmask (i -1,k,j)* tmask (i,k,j)* (rhoi (i,k,j,m) &
               - rhoi (i -1,k,j,m)) &
               + tmask (i,k,j)* tmask (i +1,k,j)* (rhoi (i +1,k,j,m) &
                                - rhoi (i,k,j,m)))                     
               e (i,k -1,j,2) = dytr (j)* p25* c1e10 &
               * (tmask (i,k -1,j -1)* tmask (i,k -1,j)* (rhoi (i,k -1,j,m) &
               - rhoi (i,k -1,j -1,m)) &
               + tmask (i,k -1,j)* tmask (i,k -1,j +1)* (rhoi (i,k -1,j +1,m)&
               - rhoi (i,k -1,j,m)) &
               + tmask (i,k,j -1)* tmask (i,k,j)* (rhoi (i,k,j,m) &
               - rhoi (i,k,j -1,m)) &
               + tmask (i,k,j)* tmask (i,k,j +1)* (rhoi (i,k,j +1,m) &
                                - rhoi (i,k,j,m)))
               e (i,k -1,j,3) = dzwr (k -1)* tmask (i,k -1,j)* tmask (i,&
                               k,j)* c1e10 &
                                * (rhoi (i,k -1,j,m) - rhoi (i,k,j,m))  
            END DO
         END DO
 
!nickbegin
!       k = km
!nickend
         DO i = 2,imm
            e (i,km,j,1) = c0
            e (i,km,j,2) = c0
            e (i,km,j,3) = c0
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     compute "K3", using "slmxr" to limit vertical slope of isopycnal
!     to guard against numerical instabilities.
!-----------------------------------------------------------------------
 
#ifdef LDD97
!lhl060506
!$OMP PARALLEL DO PRIVATE (j,k,i,chkslp,SLOPEMOD,NONDIMR,ahfctr)
      DO j = 2,jmm
         DO k = 1,km
            DO i = 1,imt
               chkslp = - sqrt (e (i,k,j,1)**2+ e (i,k,j,2)**2)* slmxr
               if (e (i,k,j,3) > chkslp) e (i,k,j,3) = chkslp
!
               SLOPEMOD= sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)/abs(e(i,k,j,3)+eps)
               F1(i,j,k)=0.5D0*( 1.0D0 + tanh((0.004D0-SLOPEMOD)/0.001D0))
               NONDIMR=-ZKP(k)/(RRD1(j)*(SLOPEMOD+eps))
               IF ( NONDIMR>=1.0 ) THEN
               F2(i,j,k)=1.0D0
               ELSE         
               F2(i,j,k) = 0.5D0*( 1.0D0 + SIN(PI*(NONDIMR-0.5D0)))
               ENDIF
               ahfctr = 1.0D0/ (e (i,k,j,3)**2+ eps)*F1(i,j,k)*F2(i,j,k)*F3(j)
!yyq 080408    ahfctr = 1.0D0/ (e (i,k,j,3)**2+ eps)*F1(i,j,k)*F2(i,j,k)
               K3 (i,k,j,1) = - e (i,k,j,3)* e (i,k,j,1)* ahfctr
               K3 (i,k,j,2) = - e (i,k,j,3)* e (i,k,j,2)* ahfctr
               K3 (i,k,j,3) = (e (i,k,j,1)**2+ e (i,k,j,2)**2)* ahfctr
            END DO
         END DO   
      END DO
!lhl060506

#else

!$OMP PARALLEL DO PRIVATE (j,k,i,chkslp,ahfctr)
      DO j = 2,jmm
         DO k = 1,km
            DO i = 2,imt -1
               chkslp = - sqrt (e (i,k,j,1)**2+ e (i,k,j,2)**2)* slmxr
               if (e (i,k,j,3) > chkslp) e (i,k,j,3) = chkslp
               ahfctr = 1.0/ (e (i,k,j,3)**2+ eps)
               K3 (i,k,j,1) = - e (i,k,j,3)* e (i,k,j,1)* ahfctr
               K3 (i,k,j,2) = - e (i,k,j,3)* e (i,k,j,2)* ahfctr
               K3 (i,k,j,3) = (e (i,k,j,1)**2+ e (i,k,j,2)**2)* ahfctr
            END DO
         END DO
      END DO
 
#endif
 
!-----------------------------------------------------------------------
!     impose zonal boundary conditions at "i"=1 and "imt"
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (j)
      DO j = 2,jmm
         CALL setbcx (K3 (1,1,j,1), imt, km)
         CALL setbcx (K3 (1,1,j,2), imt, km)
         CALL setbcx (K3 (1,1,j,3), imt, km)
      END DO
 
 
      RETURN
      END SUBROUTINE K3_123
#else
      SUBROUTINE K3_123
      RETURN
      END SUBROUTINE K3_123
#endif 
