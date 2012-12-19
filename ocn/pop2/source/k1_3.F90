!  CVS: $Id: k1_3.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     ================
      SUBROUTINE K1_3
!     ================
 
!     compute "K1(,,3)" at the center of the eastern face of "T" cells
!     use "c1e10" to keep the exponents in range.
use precision_mod 
use param_mod
use pconst_mod
use isopyc_mod
      IMPLICIT NONE
 
      REAL(r8) :: c1e10,eps,p5,c0,c1,p25,fxd,fxe,fxc,fxa,fxb,chkslp,olmask,xx
!lhl060506
      REAL(r8) , dimension(IMT,JMT,KM) :: F1,F2
      REAL(r8)  :: SLOPEMOD,NONDIMR

      F1=1.0
      F2=1.0

!lhl060506
 
!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------
 
      c1e10 = 1.0d10
 
      eps = 1.0d-25
      p5 = 0.5D0
      c0 = 0.0D0
      c1 = 1.0D0
      p25 = 0.25D0
 
!yyq 0803408
!M
!yyq 0803408
!-----------------------------------------------------------------------
!     d(rho_barx_barz)/dz centered on eastern face of "t" cells
!     Note: values involving ocean surface and ocean bottom are
!           estimated afterwards using a linear extrapolation
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,k,m,i,fxd)
      DO j = 2,jmm
         DO k = 2,km -1
            m = kisrpl (k)
            fxd = c1e10* p25* dzr (k)
            DO i = 1,imm
               e (i,k,j,3) = fxd * (rhoi (i,k -1,j,m) - rhoi (i,k +1,j,m) &
                             + rhoi (i +1,k -1,j,m) - rhoi (i +1,k +1,j,m)) 
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     linearly extrapolate densities to ocean surface for calculation
!     of d(rho_barx_barz)/dz involving level 1.
 
!     REMARK: requires min(kmt(i,jrow)) = 2 cells in ocean.
!-----------------------------------------------------------------------
 
!     k   = 1
      fxd = c1e10* dzr (1)
      fxe = dzw (0) + dzw (1)
      m = kisrpl (1)
!$OMP PARALLEL DO PRIVATE (j,i,fxa,fxb,fxc)
      DO j = 2,jmm
         DO i = 1,imm
            fxa = p5* (rhoi (i,2,j,m) + rhoi (i +1,2,j,m))
            fxb = p5* (rhoi (i,1,j,m) + rhoi (i +1,1,j,m))
            fxc = dzwr (1)* (fxb * fxe- fxa * dzw (0))
            e (i,1,j,3) = fxd * (fxc - p5* (fxa + fxb))
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     linearly extrapolate densities to ocean bottom for calculation
!     of d(rho_barx_barz)/dz involving bottom level.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 2,jmm
         DO i = 1,imm
            e (i,km,j,3) = c0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (j,i,m,k,fxa,fxb,fxc,fxe)
      DO j = 2,jmm
         DO i = 1,imm
            k = min (ITNU (i,j),ITNU (i +1,j))
            IF (k /= 0) THEN
               fxe = dzw (k -1) + dzw (k)
               m = kisrpl (k)
               fxa = p5* (rhoi (i,k -1,j,m) + rhoi (i +1,k -1,j,m))
               fxb = p5* (rhoi (i,k,j,m) + rhoi (i +1,k,j,m))
               fxc = dzwr (k -1)* (fxb * fxe- fxa * dzw (k))
               e (i,k,j,3) = dzr (k)* c1e10* (p5* (fxa + fxb) - fxc)
            END IF
 
         END DO
 
      END DO
 
 
!-----------------------------------------------------------------------
!     "e(,,,1)" = d(rho)/dx centered on east face of "T" cells
!     "e(,,,2)" = d(rho_barx_bary)/dy centered on east face of "T" cells
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i,k,m)
      DO j = 2,jmm
         DO k = 1,km
            m = kisrpl (k)
            DO i = 1,imm
               e (i,k,j,1) = tmask (i,k,j)* tmask (i +1,k,j)* OTX (j) &
                             * c1e10* (rhoi (i +1,k,j,m) - rhoi (i,k,j,m))    
               e (i,k,j,2) = p25* c1e10* dytr (J)* ( &
               rhoi (i,k,j +1,m) - rhoi (i,k,j -1,m) &
                             + rhoi (i +1,k,j +1,m) - rhoi (i +1,k,j -1,m))  
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     if any one of the 4 neighboring corner grid points is a land point
!     set "e(i,k,j,2)" to zero. note that "e(i,k,j,2)" will be used
!     only in the slope check.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i,k,olmask)
      DO j = 2,jmm
         DO k = 1,km
            DO i = 1,imm
               olmask = tmask (i,k,j -1)* tmask (i,k,j +1)* tmask (i +1,&
                       k,j -1) * tmask (i +1,k,j +1)    
               if (olmask < c1) e (i,k,j,2) = c0
            END DO
         END DO
      END DO
 
!-----------------------------------------------------------------------
!     impose zonal boundary conditions at "i"=1 and "imt"
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (j)
      DO j = 2,jmm
         CALL setbcx (e (1,1,j,1), imt, km)
         CALL setbcx (e (1,1,j,2), imt, km)
         CALL setbcx (e (1,1,j,3), imt, km)
      END DO
 
 
!lhl060506
#ifdef LDD97
!$OMP PARALLEL DO PRIVATE (j,k,i,chkslp,SLOPEMOD,NONDIMR)
      DO j = 2,jmm
         DO k = 1,km
            DO i = 1,imt
               chkslp = - sqrt (e (i,k,j,1)**2+ e (i,k,j,2)**2)* slmxr
               if (e (i,k,j,3) > chkslp) e (i,k,j,3) = chkslp
!
               SLOPEMOD= sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)/abs(e(i,k,j,3)+eps)
!
               F1(i,j,k)=0.5D0*( 1.0D0 + tanh((0.004D0-SLOPEMOD)/0.001D0))
               NONDIMR=-ZKT(k)/(RRD1(j)*(SLOPEMOD+eps))
               IF ( NONDIMR>=1.0 ) THEN
               F2(i,j,k)=1.0D0
               ELSE
               F2(i,j,k) = 0.5D0*( 1.0D0 + SIN(PI*(NONDIMR-0.5D0)))
               ENDIF
               K1 (i,k,j,3) = ( - e (i,k,j,1)* e (i,k,j,3)* &
!yyq 080408                       F1(i,j,k)*F2(i,j,k)) &
                                  F1(i,j,k)*F2(i,j,k)*F3(j)) &
                              / (e (i,k,j,3)**2+ eps)
            END DO
         END DO
      END DO

#else

!$OMP PARALLEL DO PRIVATE (j,k,i,chkslp)
      DO j = 2,jmm
         DO k = 1,km
            DO i = 1,imt
               chkslp = - sqrt (e (i,k,j,1)**2+ e (i,k,j,2)**2)* slmxr
               if (e (i,k,j,3) > chkslp) e (i,k,j,3) = chkslp
               K1 (i,k,j,3) = ( - e (i,k,j,1)* e (i,k,j,3)* fzisop (k)) &
                              / (e (i,k,j,3)**2+ eps)  
            END DO
         END DO
      END DO

#endif
!lhl060506

      RETURN
      END SUBROUTINE K1_3
 
 
#else
      SUBROUTINE K1_3
      RETURN
      END SUBROUTINE K1_3
#endif 
