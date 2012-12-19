!  CVS: $Id: convadj.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ==================
      SUBROUTINE CONVADJ
!     ==================
!     GFDL's full Convective adjustment
 
!     ---------------------------------------------------------------
!     kcon = maximum number of levels at this location
!     lcon = counts levels down
!     lcona = upper layer of a convective part of water column
!     lconb = lower layer of a convective part of water column
!     rhoup = density referenced to same level
!     rholo = density referenced to level below
!                (note that densities are not absolute!)
!     dztsum = sum of layer thicknesses
!     trasum = sum of layer tracer values
!     tramix = mixed tracer value after convection
!     lctot = total of number of levels involved in convection
!     lcven = number of levels ventilated (convection to surface)
!     note: lctot can in rare cases count some levels twice, if they
!           get involved in two originally separate, but then
!           overlapping convection areas in the water column! It
!           is a control parameter; the sensible parameter to plot
!           is lcven. Lcven is 0 on land, 1 on ocean points with no
!           convection, and anything up to km on convecting points.
!     ---------------------------------------------------------------
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use output_mod, only: ICMON
      IMPLICIT NONE
      REAL(r8)    :: RHOUP (KM),RHOLO (KM),TRASUM (2)
      REAL(r8)    :: TUP,SUP,TLO,SLO,DZTSUM,TRAMIX
      INTEGER :: KCON,LCTOT,LCVEN,L1,L,LCON,LCONA,LCONB,LMIX,n2
      REAL(r8)    :: DENS
      EXTERNAL DENS
 
      DO K = 1,KM
         RHOUP (K)= 0.0D0
         RHOLO (K)= 0.0D0
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (J,I,KCON,LCTOT,LCVEN,LCON,RHOUP,RHOLO, &
!$OMP              L,L1,TUP,SUP,TLO,SLO,K,LCONA,LCONB, &
!$OMP              DZTSUM,N,TRASUM,TRAMIX,LMIX)
      JJJ : DO J = JST,JMT
         III : DO I = 2,IMM
 
            KCON = ITNU (I,J)
            LCTOT = 0
            LCVEN = 0
            IF (KCON == 0) CYCLE III
            LCVEN = 1
            LCON = 0
 
!       FIND DENSITY OF ENTIRE ROW FOR STABILITY DETERMINATION
 
            DO L = 1,KM -1
               L1 = L +1
               TUP = AT (I,J,L1,1) - TO (L1)
               SUP = AT (I,J,L1,2) - SO (L1)
               TLO = AT (I,J, L,1) - TO (L1)
               SLO = AT (I,J, L,2) - SO (L1)
               RHOUP (L1) = DENS (TUP, SUP, L1)
               RHOLO (L) = DENS (TLO,SLO,L1)
            END DO
 
!         1. INITIAL SEARCH FOR UPPERMOST UNSTABLE PAIR; IF NONE IS
!            FOUND, MOVE ON TO NEXT COLUMN
 
            DO K = KCON -1,1, -1
               IF (RHOLO (K) > RHOUP (K +1)) LCON = K
            END DO
 
            IF (LCON == 0) CYCLE III
 
 CONV_1 : DO
            LCONA = LCON
            LCONB = LCON + 1
 
!         2. MIX THE FIRST TWO UNSTABLE LAYERS
 
            DZTSUM = DZP (LCONA) + DZP (LCONB)
            DO N = 1,2
               TRASUM (N) = AT (I,J,LCONA,N)* DZP (LCONA) + &
                            AT (I,J,LCONB,N)* DZP (LCONB) 
               TRAMIX = TRASUM (N) / DZTSUM
               AT (I,J,LCONA,N) = TRAMIX
               AT (I,J,LCONB,N) = TRAMIX
            END DO
 
!         3. TEST LAYER BELOW LCONB
 
! 1306 CONTINUE
 CONV_2 : DO

            IF (LCONB /= KCON)  THEN
            L1 = LCONB + 1
            RHOLO (LCONB) = DENS (AT (I,J,LCONB,1) - TO (L1), &
                            AT (I,J,LCONB,2) - SO (L1), L1) 
            IF (RHOLO (LCONB) > RHOUP (L1)) THEN
               LCONB = LCONB +1
               DZTSUM = DZTSUM + DZP (LCONB)
               DO N = 1,2
                  TRASUM (N) = TRASUM (N) + AT (I,J,LCONB,N)* DZP (LCONB)
                  TRAMIX = TRASUM (N) / DZTSUM
                  DO LMIX = LCONA,LCONB
                     AT (I,J,LMIX,N) = TRAMIX
                  END DO
               END DO
               CYCLE CONV_2
            END IF
            END IF
 
!         4. TEST LAYER ABOVE LCONA
 

            IF (LCONA > 1) THEN
               L1 = LCONA -1
               RHOLO (L1) = DENS (AT (I,J,L1,1) - TO (LCONA), AT (I,J,  &
                           L1,2) - SO (LCONA) ,LCONA)
               RHOUP (LCONA) = DENS (AT (I,J,LCONA,1) - TO (LCONA), AT (&
                              I,J,LCONA,2) - SO (LCONA),LCONA)
               IF (RHOLO (LCONA -1) > RHOUP (LCONA)) THEN
                  LCONA = LCONA -1
                  DZTSUM = DZTSUM + DZP (LCONA)
                  DO N = 1,2
                     TRASUM (N) = TRASUM (N) + AT (I,J,LCONA,N)* DZP (LCONA)
                     TRAMIX = TRASUM (N) / DZTSUM
                     DO LMIX = LCONA,LCONB
                        AT (I,J,LMIX,N) = TRAMIX
                     END DO
                  END DO
                  CYCLE CONV_2
               END IF
            END IF
  EXIT CONV_2
  END DO CONV_2
 
 
!         5. REMEMBER THE TOTAL NUMBER OF LEVELS MIXED BY CONVECTION
!            IN THIS WATER COLUMN, AS WELL AS THE VENTILATED COLUMN
 
            LCTOT = LCTOT + LCONB - LCONA + 1
            IF (LCONA == 1) LCVEN = LCONB - LCONA + 1
 
!         6. RESUME SEARCH IF STEP 3. AND 4. HAVE BEEN PASSED AND THIS
!            UNSTABLE PART OF THE WATER COLUMN HAS THUS BEEN REMOVED,
!            I.E. FIND FURTHER UNSTABLE AREAS FURTHER DOWN THE COLUMN
 
            IF (LCONB == KCON) THEN
            ICMON (I,J,1)= ICMON (I,J,1) + LCTOT
            ICMON (I,J,2)= ICMON (I,J,2) + LCVEN
            CYCLE III
            ENDIF
            LCON = LCONB
 
 CONV_3 : DO

            LCON = LCON + 1
            IF (LCON == KCON) THEN 
            ICMON (I,J,1)= ICMON (I,J,1) + LCTOT
            ICMON (I,J,2)= ICMON (I,J,2) + LCVEN
            CYCLE III
            ENDIF
            IF (RHOLO (LCON) <= RHOUP (LCON +1)) CYCLE CONV_3
            EXIT CONV_3
  END DO CONV_3
!

  END DO CONV_1
 
 
!JXZ------------------------------------------------------------------
!            ICMON (I,J,1)= ICMON (I,J,1) + LCTOT
!            ICMON (I,J,2)= ICMON (I,J,2) + LCVEN
!JXZ------------------------------------------------------------------
 
         END DO III
      END DO JJJ
     
      if (nx_proc.ne.1) then
      do n =1,2
      call exchange_2d_real(icmon(1,1,n),1,0)
      end do
      n2=mod(mytid,nx_proc)
      do j=1,jmt
        if (n2==0) then
         icmon(1,j,1)=0.0
         icmon(1,j,2)=0.0
        end if
        if (n2==nx_proc-1)then
         icmon(imt,j,1)=0.0 
         icmon(imt,j,2)=0.0 
        end if
      end do
      end if
! 
!JXZ------------------------------------------------------------------
!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!JXZ------------------------------------------------------------------
 
      DO N = 1,2
!$OMP PARALLEL DO PRIVATE (K,J)
         DO K = 1,KM
            DO J = JST,JMT
               AT (1,J,K,N) = AT (IMM,J,K,N)
               AT (IMT,J,K,N) = AT (2,J,K,N)
            END DO
         END DO
      call exch_boundary(at(1,1,1,n),km)
      END DO
      
#ifdef TSPAS
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = JST,JMT
            DO I = 1  ,IMT
                ATB(I,J,K,1)  =  AT(I,J,K,1)
                ATB(I,J,K,2)  =  AT(I,J,K,2)
            END DO
            END DO
         call exch_boundary(ATB(1,1,1,1),km)
         call exch_boundary(ATB(1,1,1,2),km)
         END DO
#endif
      RETURN
      END SUBROUTINE CONVADJ
