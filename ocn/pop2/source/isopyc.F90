!  CVS: $Id: isopyc.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     =================
      SUBROUTINE ISOPYC
!     =================
 
!     Compute the isopycnal mixing tensor components and the
!     isopycnal advection velocities which param_modeterize the effect
!     of eddies on the isopycnals.
 
!     Mixing tensor "K" is ...
!          | 1.0            K1(,,,2)          K1(,,,3) |
!          |                                           |
!     K =  | K2(,,,1)        1.0              K2(,,,3) |
!          |                                           |
!          | K3(,,,1)       K3(,,,2)          K3(,,,3) |
 
!     where K1(,,,2) and K2(,,,1) are set to 0.0 (neglected)
 
!     output:
!       rhoi = density at tau-1 referenced to pressure levels
!       K1   = tensor components (1,2), (1,3) centered on east face
!              of "T" cells
!       K2   = tensor components (2,1), (2,3) centered on north face
!              of "T" cells
!       K3   = tensor components (3,1), (3,2), (3,3) centered on
!              bottom face of "T" cells
!       adv_vetiso = isopycnal advective vel on east face of "T" cell
!       adv_vntiso = isopycnal advective vel on north face of "T" cell
!               (Note: this includes the cosine factor as in "adv_vnt")
!       adv_vbtiso = isopycnal advective vel on bottom face of "T" cell
use precision_mod 
use param_mod
use pconst_mod
use tracer_mod
use isopyc_mod
      IMPLICIT NONE
 
      REAL(r8):: tq,sq,tref0,sref0
      INTEGER :: kq
      REAL(r8):: DENS
      EXTERNAL DENS
 
!-----------------------------------------------------------------------
!     compute normalized densities for each isopycnal reference pressure
!     level using a 3rd order polynomial fit to the equation of state.
!     for each isopycnal reference pressure level, the same reference
!     potential temperature, reference salinity and expansion coeff
!     values are used at all of the vertical levels.
 
!     Note: this density is used for the mixing tensor in both the
!     Redi/Cox and Gent/McWilliams options
!-----------------------------------------------------------------------
 
      allocate(e(imt,kmp1,jmt,3),rhoi(imt,0:km,jmt,nrpl))
      allocate(K1(imt,0:km,jmt,3:3),K2(imt,0:km,jmt,3:3),K3(imt,0:km,jmt,1:3))


!$OMP PARALLEL DO PRIVATE (j,k,i,m)
      DO m = 1,nrpl
         DO j = 1,jmt
!           DO k = 0,km
               DO i = 1,imt
                  rhoi (i,0,j,m) = 0.0D0
               END DO
!           END DO
         END DO
      END DO


      DO m = 1,nrpl
         tref0 = to (krplin (m))
         sref0 = so (krplin (m))
!$OMP PARALLEL DO PRIVATE (j,k,i,tq,sq,kq)
       DO k = 1,km
          DO j = 1,jmt
               DO i = 1,imt
                  tq = atb (i,j,k,1) - tref0
                  sq = atb (i,j,k,2) - sref0
                  kq = krplin (m)
                  rhoi (i,k,j,m) = dens (tq, sq, kq)
               END DO
            END DO
         END DO
      END DO

 
!-----------------------------------------------------------------------
!     evaluate K2(,,3) centered on the northern face of "T" cells
!-----------------------------------------------------------------------
 
      CALL k2_3
 
!-----------------------------------------------------------------------
!     evaluate K1(,,3) centered on eastern face of "T" cells
!-----------------------------------------------------------------------
 
      CALL k1_3
 
!-----------------------------------------------------------------------
!     evaluate K3(,,1..3) centered on bottom face of "T" cells
!-----------------------------------------------------------------------
 
      CALL k3_123
 
!-----------------------------------------------------------------------
!     compute isopycnal advective velocities for tracers
!-----------------------------------------------------------------------
 
      deallocate(e,rhoi)
      allocate(adv_vetiso(imt,km,jmt),adv_vntiso(imt,km,jmt),adv_vbtiso(imt,0:km,jmt))

      CALL isoadv
 
      RETURN
      END SUBROUTINE ISOPYC
 
 
#else
      SUBROUTINE ISOPYC ()
      RETURN
 
      END SUBROUTINE ISOPYC
 
#endif 
