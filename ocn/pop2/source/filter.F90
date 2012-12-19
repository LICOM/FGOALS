!  CVS: $Id: filter.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =======================
      SUBROUTINE FILTER (X,Z,KK,JFN1,JFN2)
!     =======================
!     1-D zonal fouriour filtering.
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
      IMPLICIT NONE
      INTEGER :: JFS1,JFS2,JFN1,JFN2,KK, NNS, NNN, IMM2
      REAL(r8)    :: XM,A1,A2,B1,B2,OMM
      REAL(r8)    :: X (IMT,JMT,KK),XS (IMT),Z (IMT,JMT,KM)
 
      OMM = 1.0D0/ FLOAT (IMT -2)
 
!$OMP PARALLEL DO PRIVATE (K,J,I,XS,XM,A1,B1,A2,B2)
      DO K = 1,KK
 
         DO J = JFN1,JFN2
 
            DO I = 1,IMT
               XS (I)= X (I,J,K)* Z (I,J,K)
            END DO
 
 
            XM = 0.0D0
            A1 = 0.0D0
            B1 = 0.0D0
            A2 = 0.0D0
            B2 = 0.0D0
 
            DO I = 2,IMM
               XM = XM + XS (I)* OMM
               A1 = A1+ XS (I)* CF1 (I)* OMM
               B1 = B1+ XS (I)* SF1 (I)* OMM
               A2 = A2+ XS (I)* CF2 (I)* OMM
               B2 = B2+ XS (I)* SF2 (I)* OMM
            END DO
 
 
            DO I = 2,IMM
               XS (I)= XM + A1* CF1 (I) + B1* SF1 (I) + A2* CF2 (I)     &
                    + B2* SF2 (I)
               X (I,J,K)= XS (I)* Z (I,J,K)
            END DO
 
            X (1,J,K)= X (IMM,J,K)
            X (IMT,J,K)= X (2,J,K)
 
!YU   NNN=2*(7-J)
 
!YU   DO 102 N=1,NNN
!YU
!YU   DO I=1,IMT
!YU   XS(I)=X(I,J,K)*Z(I,J,K)
!YU   ENDDO
 
!YU   DO I=2,IMM
!YU   X(I,J,K)=(0.5*XS(I)+0.25*(XS(I-1)+XS(I+1)))*Z(I,J,K)
!YU   ENDDO
 
!YU   X(1  ,J,K)=X(IMM,J,K)
!YU   X(IMT,J,K)=X(2  ,J,K)
!YU
!YU  102 CONTINUE
!YU
         END DO
 
      END DO
 
      RETURN
      END SUBROUTINE FILTER
 
 
