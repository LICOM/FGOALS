!  CVS: $Id: density.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ==================
      SUBROUTINE DENSITY
!     ==================
 
!     computes normalized densities by using a 3rd order polynomial
!     fit to the equation of state set by Joint Panel on Oceanographic
!     Tables & Standards (UNESCO, 1981).
 
!     note.. for precision purposes, there is a depth dependent
!     constant subtracted from the density returned by this routine.
!     so... this routine should be used only for horizontal gradients
!     of density.
 
!     inputs:
 
!     t  = temperatures (potential deg C)
!     s  = salinities   (ppt-35)/1000)
 
!     output:
 
!     rho = normalized densities
 
!     These densities represent the in situ density at a level minus
!     a depth dependent normalization. The complete in situ density is
!     rho_complete = dens(t(k)-to(k),s(k)-so(k),k) + rho_norm(k)
!     where rho_norm(k) are the depth dependent normalization densities
!     given at the bottom of dncoef.h
 
!     dRHO=c1*t+c2*s+c3*tt+c4*ts*c5*ss+c6*ttt+c7*tss+c8*tts+c9*sss

use precision_mod 
use param_mod
use pconst_mod
use tracer_mod
use work_mod, only: wka
      IMPLICIT NONE
 
 
      REAL(r8)    :: TQ,SQ
 
!$OMP PARALLEL DO PRIVATE (K,J,I,TQ,SQ) 
      DO K = 1,KM
!YU   DO J=1,JMT
!lhl      DO J=JST,JMT
         DO J = JST,JET
            DO I = 1,IMT
               IF (VIT (I,J,K) > 0.0) THEN
                  TQ = AT (I,J,K,1) - TO (K)
                  SQ = AT (I,J,K,2) - SO (K)
!lhl0608                  TQ = ATB(I,J,K,1) - TO (K)
!lhl0608                  SQ = ATB(I,J,K,2) - SO (K)
!lhl1204
                  PDENSITY(I,J,K)=1.0d+3+PO(K)+(C(K,1)+(C(K,4)+C(K,7)*SQ)*SQ&
                                   +(C(K,3)+C(K,8)*SQ+C(K,6)*TQ)*TQ)*TQ&
                                   +(C(K,2)+(C(K,5)+C(K,9)*SQ)*SQ)*SQ
!                  WKA(I,J,K)=      (C(K,1)+(C(K,4)+C(K,7)*SQ)*SQ+(C(K,3)+C(K,8)*SQ+&
!                  C(K,6)*TQ)*TQ)*TQ+(C(K,2)+(C(K,5)+C(K,9)*SQ)*SQ)*SQ
!lhl1204
               ELSE
!lhl1204
!                  WKA (I,J,K) = 0.0
                  PDENSITY (I,J,K) = 0.0
!lhl1204
               END IF
            END DO
         END DO
      END DO
      RETURN
      END SUBROUTINE DENSITY
 
 
