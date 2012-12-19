!  CVS: $Id: smth.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =========================
      SUBROUTINE SMOOTH (X,IT,S)
!     =========================
!   9-point smooth
use precision_mod
use param_mod
      IMPLICIT NONE
      REAL(r8)    :: X (imt,jmt),XS (imt,jmt)
      INTEGER :: it (imt,jmt)
      INTEGER :: IT5,IT6,IT7,IT8,IT1,IT2,IT3,IT4
      REAL(r8)    :: XS5,XS6,XS7,XS8,XS1,XS2,XS3,XS4
      REAL(r8)    :: S1,BT0,BT1,BT2,S
      S1 = 1.0D0- S
 
      BT0 = S1* S1
      BT1 = 0.5D0* S * S1
      BT2 = 0.25D0* S * S
      DO J = 1,jmt
         DO I = 1,imt
            XS (I,J)= X (I,J)
         END DO
 
      END DO
 
      JJJ : DO J = 2,JMM
         III : DO I = 2,IMM
            IF (IT (I,J) == 0) CYCLE III
            IT5 = IT (I +1,J +1)
            XS5 = XS (I +1,J +1)
            IT6 = IT (I -1,J +1)
            XS6 = XS (I -1,J +1)
            IT7 = IT (I -1,J -1)
            XS7 = XS (I -1,J -1)
            IT8 = IT (I +1,J -1)
            XS8 = XS (I +1,J -1)
            IT1 = IT (I +1,J)
            XS1 = XS (I +1,J)
            IT2 = IT (I -1,J)
            XS2 = XS (I -1,J)
            IT3 = IT (I,J +1)
            XS3 = XS (I,J +1)
            IT4 = IT (I,J -1)
            XS4 = XS (I,J -1)
            IF (IT1 == 0) XS1 = XS (I,J)
            IF (IT2 == 0) XS2 = XS (I,J)
            IF (IT3 == 0) XS3 = XS (I,J)
            IF (IT4 == 0) XS4 = XS (I,J)
            IF (IT5 == 0) XS5 = XS (I,J)
            IF (IT6 == 0) XS6 = XS (I,J)
            IF (IT7 == 0) XS7 = XS (I,J)
            IF (IT8 == 0) XS8 = XS (I,J)
            X (I,J)= BT0* XS (I,J) + BT1* (XS1+ XS2+ XS3+ XS4) + BT2* ( &
                 XS5+ XS6+ XS7+ &
            XS8)
         END DO III
      END DO JJJ
      DO J = 2,JMM
         X (1,J) = X (IMM,J)
         X (IMT,J) = X (2,J)
      END DO
 
      RETURN
      END SUBROUTINE SMOOTH
 
 
 
!     =========================
      SUBROUTINE SMOOTH2 (X,IT,S,FS1,FE1,FS2,FE2)
!     =========================
!   9-point smooth
use precision_mod
use param_mod
      IMPLICIT NONE
      REAL(r8)    :: X (imt,jmt),XS (imt,jmt)
      INTEGER :: it (imt,jmt),FS (2),FE (2),NN
      INTEGER :: FS1,FS2,FE1,FE2
      INTEGER :: IT5,IT6,IT7,IT8,IT1,IT2,IT3,IT4
      REAL(r8)    :: XS5,XS6,XS7,XS8,XS1,XS2,XS3,XS4
      REAL(r8)    :: S1,BT0,BT1,BT2,S
      FS (1)= FS1
 
      FE (1)= FE1
      FS (2)= FS2
      FE (2)= FE2
      S1 = 1.0- S
      BT0 = S1* S1
      BT1 = 0.5D0* S * S1
      BT2 = 0.25D0* S * S
      DO J = 1,jmt
         DO I = 1,imt
            XS (I,J)= X (I,J)
         END DO
 
      END DO
 
      DO NN = 1,2
         JJJ : DO J = FS (NN),FE (NN)
            III : DO I = 2,IMM
               IF (IT (I,J) == 0) CYCLE III
               IT5 = IT (I +1,J +1)
               XS5 = XS (I +1,J +1)
               IT6 = IT (I -1,J +1)
               XS6 = XS (I -1,J +1)
               IT7 = IT (I -1,J -1)
               XS7 = XS (I -1,J -1)
               IT8 = IT (I +1,J -1)
               XS8 = XS (I +1,J -1)
               IT1 = IT (I +1,J)
               XS1 = XS (I +1,J)
               IT2 = IT (I -1,J)
               XS2 = XS (I -1,J)
               IT3 = IT (I,J +1)
               XS3 = XS (I,J +1)
               IT4 = IT (I,J -1)
               XS4 = XS (I,J -1)
               IF (IT1 == 0) XS1 = XS (I,J)
               IF (IT2 == 0) XS2 = XS (I,J)
               IF (IT3 == 0) XS3 = XS (I,J)
               IF (IT4 == 0) XS4 = XS (I,J)
               IF (IT5 == 0) XS5 = XS (I,J)
               IF (IT6 == 0) XS6 = XS (I,J)
               IF (IT7 == 0) XS7 = XS (I,J)
               IF (IT8 == 0) XS8 = XS (I,J)
               X (I,J)= BT0* XS (I,J) + BT1* (XS1+ XS2+ XS3+ XS4)       &
                     + BT2* (XS5+ XS6+ &
               XS7+ XS8)
            END DO III
         END DO JJJ
         DO J = 2,JMM
            X (1,J) = X (IMM,J)
            X (IMT,J) = X (2,J)
         END DO
 
      END DO
      RETURN
      END SUBROUTINE SMOOTH2
 
 
 
!     =========================
      SUBROUTINE SMOOTH3 (X,VT,S,FS1,FE1,FS2,FE2)
!     =========================
!   9-point smooth
use param_mod
      IMPLICIT NONE
      REAL    :: X (imt,jmt,km),VT (imt,jmt,km),XS (imt,jmt)
      INTEGER :: it (imt,jmt),FE (2),FS (2),NN
      INTEGER :: FS1,FS2,FE1,FE2
      INTEGER :: IT5,IT6,IT7,IT8,IT1,IT2,IT3,IT4
      REAL    :: XS5,XS6,XS7,XS8,XS1,XS2,XS3,XS4
      REAL    :: S1,BT0,BT1,BT2,S
      FS (1)= FS1
 
      FE (1)= FE1
      FS (2)= FS2
      FE (2)= FE2
      S1 = 1.0D0- S
      BT0 = S1* S1
      BT1 = 0.5D0* S * S1
      BT2 = 0.25D0* S * S
 
      DO K = 1,KM
         DO J = 1,jmt
            DO I = 1,imt
               IT (I,J)= int (VT (I,J,K))
            END DO
 
         END DO
 
         DO J = 1,jmt
            DO I = 1,imt
               XS (I,J)= X (I,J,K)
            END DO
 
         END DO
 
         DO NN = 1,2
            JJJ : DO J = FS (NN),FE (NN)
               III : DO I = 2,IMM
                  IF (IT (I,J) == 0) CYCLE III
                  IT5 = IT (I +1,J +1)
                  XS5 = XS (I +1,J +1)
                  IT6 = IT (I -1,J +1)
                  XS6 = XS (I -1,J +1)
                  IT7 = IT (I -1,J -1)
                  XS7 = XS (I -1,J -1)
                  IT8 = IT (I +1,J -1)
                  XS8 = XS (I +1,J -1)
                  IT1 = IT (I +1,J)
                  XS1 = XS (I +1,J)
                  IT2 = IT (I -1,J)
                  XS2 = XS (I -1,J)
                  IT3 = IT (I,J +1)
                  XS3 = XS (I,J +1)
                  IT4 = IT (I,J -1)
                  XS4 = XS (I,J -1)
                  IF (IT1 == 0) XS1 = XS (I,J)
                  IF (IT2 == 0) XS2 = XS (I,J)
                  IF (IT3 == 0) XS3 = XS (I,J)
                  IF (IT4 == 0) XS4 = XS (I,J)
                  IF (IT5 == 0) XS5 = XS (I,J)
                  IF (IT6 == 0) XS6 = XS (I,J)
                  IF (IT7 == 0) XS7 = XS (I,J)
                  IF (IT8 == 0) XS8 = XS (I,J)
                  X (I,J,K)= BT0* XS (I,J) + BT1* (XS1+ XS2+ XS3+ XS4)  &
                        + BT2* (XS5+ XS6+ XS7+ XS8)
               END DO III
            END DO JJJ
            DO J = 2,JMM
               X (1,J,K) = X (IMM,J,K)
               X (IMT,J,K) = X (2,J,K)
            END DO
 
         END DO
      END DO
 
      RETURN
      END SUBROUTINE SMOOTH3
 
 
