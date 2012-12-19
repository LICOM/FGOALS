!  CVS: $Id: bclinc.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE BCLINC
!     =================
!     INTEGRATION OF MOMENTUM EQUATION
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use forc_mod, only: psa,su,sv
      IMPLICIT NONE
      REAL(r8)    :: AA,GGU,WK1,WK2,fil_lat1,fil_lat2

!lhl0711
#ifdef CANUTO
      REAL(r8)    :: AIDIF
      REAL(r8)    :: SBCX(imt,jmt),BBCX(imt,jmt)
      REAL(r8)    :: SBCY(imt,jmt),BBCY(imt,jmt)


!M
!$OMP PARALLEL DO PRIVATE (J,I) 
         DO J = JSM,JEM
            DO I = 2,IMM
               SBCX(I,J) = SU (I,J)* OD0
               SBCY(I,J) = SV (I,J)* OD0
               BBCX(I,J)= C0F*SQRT(UP(I,J,KM)*UP(I,J,KM)+VP(I,J,KM)*VP(I,J,KM))&
                          *(UP(I,J,KM)*CAG+SNLAT(J)*VP (I,J,KM)*SAG)
               BBCY(I,J)= C0F*SQRT(UP(I,J,KM)*UP(I,J,KM)+VP(I,J,KM)*VP(I,J,KM))&
                          *(-SNLAT(J)*UP(I,J,KM)*SAG+VP(I,J,KM)*CAG)
            ENDDO
         ENDDO

      AIDIF=0.0  
!      if (ISC/=0)  AIDIF=0.5
 
#endif
!
!---------------------------------------------------------------------
!      Define the threthold latitute for zonal smoother
       fil_lat1=63.0D0
       fil_lat2=63.0D0
!lhl0711
!---------------------------------------------------------------------
!     ADVECTION + DIFFUSION + CORIOLIS
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               DLU (I,J,K) = DLU (I,J,K) - FF (J)* VP (I,J,K)
               DLV (I,J,K) = DLV (I,J,K) + FF (J)* UP (I,J,K)
            END DO
         END DO
      END DO
 
!---------------------------------------------------------------------
!     PRESSURE GRADIENT FORCES
!---------------------------------------------------------------------
 
!!@@@@ DP'/DX
 
      AA = 0.0
      IF (ISC /= 0) AA = 0.5
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            H0BF (I,J)= H0BF (I,J)* ONBB
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            WORK (I,J)= AA * H0BF (I,J) + (1.0- AA)* H0BL (I,J)
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (J,I,WKK)
      DO J = JST,JET
         DO I = 1,IMT
            WKK (1)= (PSA (I,J)* OD0+ WORK (I,J)* G)* VIT (I,J,1)
            DO K = 1,KM
               WKK (K +1)= WKK (K) - GG (I,J,K)* DZP (K)* VIT (I,J,K)
            END DO
 
            DO K = 1,KM
               WKA (I,J,K)= 0.25* (WKK (K) + WKK (K +1))
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               DLV (I,J,K) = DLV (I,J,K) - OUY (J)* &
                             (WKA (I,J +1,K) - WKA (I,J,K) + WKA (I -1, &
                              J +1,K) - WKA (I -1,J,K))    
               DLU (I,J,K) = DLU (I,J,K) - OUX (J)* &
                             (WKA (I,J,K) - WKA (I -1,J,K) + WKA (I,    &
                              J +1,K) - WKA (I -1,J +1,K))
            END DO
         END DO
      END DO
 
!!@@@@ G'DH/DX
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,K)= (1.0+ OHBT (I,J)* ZKT (K))* WORK (I,J)      &
                            * VIT (I,J,K)
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I,GGU)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               GGU = 0.125* (GG (I,J,K) + GG (I -1,J,K) + GG (I,J +1,K) &
                     + GG (I -1,J +1,K))
               DLV (I,J,K) = DLV (I,J,K) + GGU * OUY (J)* VIV (I,J,K)* &
                             (WKA (I,J +1,K) - WKA (I,J,K) + WKA (I -1, &
                              J +1,K) - WKA (I -1,J,K)) 
               DLU (I,J,K) = DLU (I,J,K) + GGU * OUX (J)* VIV (I,J,K)* &
                             (WKA (I,J,K) - WKA (I -1,J,K) + WKA (I,    &
                              J +1,K) - WKA (I -1,J +1,K)) 
            END DO
         END DO
      END DO

 
!---------------------------------------------------------------------
!     CORIOLIS ADJUSTMENT
!---------------------------------------------------------------------
 
      IF (ISC == 0) THEN
!$OMP PARALLEL DO PRIVATE (K,J,I,WK1,WK2)
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 2,IMM
                  WK1 = EPEA (J)* DLV (I,J,K) + EPEB (J)* DLU (I,J,K)
                  WK2 = EPEA (J)* DLU (I,J,K) - EPEB (J)* DLV (I,J,K)
                  DLV (I,J,K)= WK1* VIV (I,J,K)
                  DLU (I,J,K)= WK2* VIV (I,J,K)
               END DO
            END DO
         END DO
      ELSE
!$OMP PARALLEL DO PRIVATE (K,J,I,WK1,WK2)
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 2,IMM
                  WK1 = EPLA (J)* DLV (I,J,K) + EPLB (J)* DLU (I,J,K)
                  WK2 = EPLA (J)* DLU (I,J,K) - EPLB (J)* DLV (I,J,K)
                  DLV (I,J,K)= WK1* VIV (I,J,K)
                  DLU (I,J,K)= WK2* VIV (I,J,K)
               END DO
            END DO
         END DO
      END IF
 
 
!---------------------------------------------------------------------
!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!---------------------------------------------------------------------
 
      if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1,KM
         DO J = JSM,JEM
            DLV (1,J,K) = DLV (IMM,J,K)
            DLU (1,J,K) = DLU (IMM,J,K)
            DLV (IMT,J,K) = DLV (2,J,K)
            DLU (IMT,J,K) = DLU (2,J,K)
         END DO
      END DO
      end if
!
 
#ifdef SPMD
      call exch_boundary(dlu,km)
      call exch_boundary(dlv,km)
#endif
!YU  Oct. 24, 2005

      CALL SMUV (DLU,VIV,KM,fil_lat1)
      CALL SMUV (DLV,VIV,KM,fil_lat1)
!YU 

 
!---------------------------------------------------------------------
!     PREDICTING VC & UC
!---------------------------------------------------------------------
 
      IF (ISC < 1) THEN
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               V (I,J,K)= VP (I,J,K) + DLV (I,J,K)* DTC
               U (I,J,K)= UP (I,J,K) + DLU (I,J,K)* DTC
            END DO
         END DO
      END DO
 
!---------------------------------------------------------------------
!@@@  INTERACTION BETWEEN BAROTROPIC AND BAROCLINIC MODES
!---------------------------------------------------------------------
      CALL VINTEG (U,WORK)
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               U (I,J,K)= (U (I,J,K) - WORK (I,J) + UB (I,J))* VIV (I,J,K)
            END DO
         END DO
      END DO
!
#ifdef SPMD
      call exch_boundary(u,km)
#endif
 
      CALL VINTEG (V,WORK)
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               V (I,J,K)= (V (I,J,K) - WORK (I,J) + VB (I,J))* VIV (I,J,K)
            END DO
         END DO
      END DO

#ifdef SPMD
      call exch_boundary(v,km)
#endif
 
 
      ISC = ISC +1
    
      ELSE 

 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,K)= VP (I,J,K) + DLV (I,J,K)* DTC2
            END DO
         END DO
      END DO

!lhl0711
!#ifdef CANUTO
!      CALL INVTRI1 (WKA,SBCY,BBCY,AKMU,AIDIF,DTC2)
!#endif
!lhl0711
 
      CALL VINTEG (WKA,WORK)
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,K)= (WKA (I,J,K) - WORK (I,J) + VB (I,J))* VIV (I,J,K)
            END DO
         END DO
      END DO
 
!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------
 
#ifdef SPMD
      call exch_boundary(wka,km)
#endif

 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               VP (I,J,K) = AFC2* V (I,J,K) + AFC1* (VP (I,J,K) + WKA (I,J,K))
               V (I,J,K) = WKA (I,J,K)
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,K)= UP (I,J,K) + DLU (I,J,K)* DTC2
            END DO
         END DO
      END DO
 
!lhl0711
!#ifdef CANUTO
!      CALL INVTRI1 (WKA,SBCX,BBCX,AKMU,AIDIF,DTC2)
!#endif
!lhl0711
 
      CALL VINTEG (WKA,WORK)
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,K)= (WKA (I,J,K) - WORK (I,J) + UB (I,J))* VIV (I,J,K)
            END DO
         END DO
      END DO
 
!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------
 
#ifdef SPMD
      call exch_boundary(wka,km)
#endif
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UP (I,J,K) = AFC2* U (I,J,K) + AFC1* (UP (I,J,K) + WKA (I,J,K))
               U (I,J,K) = WKA (I,J,K)
            END DO
         END DO
      END DO
 
!YU  Oct. 24, 2005
!lhl0711     IF (MOD(ISC,60)==0) THEN
     IF (MOD(ISC,160)==1) THEN
        CALL SMUV (U,VIV,KM,fil_lat2)
        CALL SMUV (V,VIV,KM,fil_lat2)
        CALL SMUV (UP,VIV,KM,fil_lat2)
        CALL SMUV (VP,VIV,KM,fil_lat2)
#ifdef SPMD
        call exch_boundary(u,km)
        call exch_boundary(v,km)
        call exch_boundary(up,km)
        call exch_boundary(vp,km)
#endif
     END IF
!YU 
 
      ISC = ISC +1
      END IF

 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UTF (I,J,K)= UTF (I,J,K) + U (I,J,K)
               VTF (I,J,K)= VTF (I,J,K) + V (I,J,K)
            END DO
         END DO
      END DO
 
!
      RETURN
      END SUBROUTINE BCLINC
