!  CVS: $Id: dens.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =======================
      FUNCTION DENS (TQ,SQ,KK)
!     =======================
use precision_mod 
use param_mod
use pconst_mod
      IMPLICIT NONE
 
      REAL(r8)    :: DENS,TQ,SQ
      integer:: KK 
      DENS = (C (KK,1) + (C (KK,4) + C (KK,7)* SQ)* SQ + &
      (C (KK,3) + C (KK,8)* SQ + C (KK,6)* TQ)* TQ)* TQ + &
             (C (KK,2) + (C (KK,5) + C (KK,9)* SQ)* SQ)* SQ  
 
      RETURN
      END FUNCTION DENS
 
 
