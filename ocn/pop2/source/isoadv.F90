!  CVS: $Id: isoadv.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     =================
      SUBROUTINE ISOADV
!     =================
 
!     compute isopycnal transport velocities.
use precision_mod 
use param_mod
use pconst_mod
use isopyc_mod
#ifdef SPMD
use msg_mod
#endif

      IMPLICIT NONE
 
      REAL(r8):: p5,c0,fxa
      INTEGER :: jstrt
 
!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------
      p5 = 0.5D0
 
      c0 = 0.0D0
 
!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal mixing velocity
!     at the center of the northern face of the "t" cells.
!-----------------------------------------------------------------------
      adv_vntiso=c0 
      adv_vbtiso=c0 
      adv_vetiso=c0 
!$OMP PARALLEL DO PRIVATE (j,k,i,fxa)
      DO j = 2,jmm
         DO k = 2,km -1
            fxa = - p5* dzr (k)* athkdf * SINU (j)
            DO i = 1,imt
               adv_vntiso (i,k,j) = fxa * tmask (i,k,j)* tmask (i,k,j +1)* ( &
                                    K2 (i,k -1,j,3) - K2 (i,k +1,j,3))
            END DO
         END DO
      END DO
 
 
!     consider the top and bottom levels. "K2" is assumed to be zero
!     at the ocean top and bottom.
 
!     k = 1
      fxa = - p5* dzr (1)* athkdf
!$OMP PARALLEL DO PRIVATE (j,i,k)
      DO j = 2,jmm
         DO i = 1,imt
            adv_vntiso (i,1,j) = - fxa * tmask (i,1,j)* tmask (i,1,j +1)&
                                * SINU (j) &
                                 * (K2 (i,1,j,3) + K2 (i,2,j,3))   
         END DO
      END DO
 
 
 
 
!$OMP PARALLEL DO PRIVATE (j,i,k)
      DO j = 2,jmm
         DO i = 1,imt
            k = min (ITNU (i,j),ITNU (i,j +1))
            IF (k /= 0) THEN
               adv_vntiso (i,k,j) = - p5* dzr (k)* athkdf * SINU (j)    &
                                   * tmask (i,k,j) &
                                    * tmask (i,k,j +1)* (K2 (i,k,j,3)   &
                                      + K2 (i,k -1,j,3))          
            END IF
         END DO
      END DO

#ifdef SPMD
      call exchange_3d_iso(adv_vntiso,km,1,1)
#endif 
 
!-----------------------------------------------------------------------
!     compute the zonal component of the isopycnal mixing velocity
!     at the center of the eastern face of the "t" grid box.
!-----------------------------------------------------------------------
 
      jstrt = 2
 
!$OMP PARALLEL DO PRIVATE (j,k,fxa,i)
      DO j = jstrt,jmm
         DO k = 2,km -1
            fxa = - p5* dzr (k)* athkdf
            DO i = 1,imm
               adv_vetiso (i,k,j) = fxa * tmask (i,k,j)* tmask (i +1,k,j) &
                                    * (K1 (i,k -1,j,3) - K1 (i,k +1,j,3))      
            END DO
         END DO
      END DO
 
 
!     consider the top and bottom levels. "K1" is assumed to be zero
!     at the ocean top and bottom.
 
!     k = 1
      fxa = - p5* dzr (1)* athkdf
!$OMP PARALLEL DO PRIVATE (j,k,i)
      DO j = jstrt,jmm
         DO i = 1,imm
            adv_vetiso (i,1,j) = - fxa * tmask (i,1,j)* tmask (i +1,1,j) &
                                 * (K1 (i,1,j,3) + K1 (i,2,j,3)) 
         END DO
      END DO
 
 
 
 
!$OMP PARALLEL DO PRIVATE (j,i,k)
      DO j = jstrt,jmm
         DO i = 1,imm
            k = min (ITNU (i,j),ITNU (i +1,j))
            IF (k /= 0) THEN
               adv_vetiso (i,k,j) = - p5* dzr (k)* athkdf * tmask (i,k,j) &
                                    * tmask (i +1,k,j)* (K1 (i,k,j,3)   &
                                      + K1 (i,k -1,j,3)) 
            END IF
         END DO
      END DO
 
 
!----------------------------------------------------------------------
!     set the boundary conditions
!----------------------------------------------------------------------
      if (nx_proc ==1) then
!$OMP PARALLEL DO PRIVATE (j)
      DO j = jstrt,jmm
         CALL setbcx (adv_vetiso (1,1,j), imt, km)
      END DO
      end if
#ifdef SPMD
      call exchange_3d_iso(adv_vntiso,km,1,1)
#endif
 
!----------------------------------------------------------------------
!     compute the vertical component of the isopycnal mixing velocity
!     at the center of the bottom face of the "t" cells, using the
!     continuity equation for the isopycnal mixing velocities
!-----------------------------------------------------------------------
 
 
 
!$OMP PARALLEL DO PRIVATE (j,i,k)
      DO j = jstrt,jmm
         DO k = 1,km -1
            DO i = 2,imt
               adv_vbtiso (i,k,j) = DZP (k)* ( &
               (adv_vetiso (i,k,j) - adv_vetiso (i -1,k,j))* OTX (j) + &
                                    (adv_vntiso (i,k,j) - adv_vntiso (i,&
                                     k,j -1))* cstrdytr (j))     
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (j)
      DO j = jstrt,jmm
         DO k = 1,km -1
            DO i = 2,imt
               adv_vbtiso (i,k,j) = adv_vbtiso (i,k,j) + adv_vbtiso (i,k -1,j)
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = jstrt,jmm
         DO i = 2,imt
            adv_vbtiso (i,ITNU (i,j),j) = c0
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     set the boundary conditions
!-----------------------------------------------------------------------
      if (nx_proc ==1) then
!$OMP PARALLEL DO PRIVATE (j)
      DO j = jstrt,jmm
         CALL setbcx (adv_vbtiso (1,0,j), imt, km +1)
      END DO
      end if 
#ifdef SPMD
      call exchange_3d_iso(adv_vntiso,km,1,0)
#endif

      RETURN
      END SUBROUTINE ISOADV
 
 
#else
      SUBROUTINE ISOADV
      RETURN
      END SUBROUTINE ISOADV
#endif 
