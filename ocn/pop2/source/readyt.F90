!  CVS: $Id: readyt.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE READYT
!     =================
!     PREPARATION OF BAROTROPIC AND BAROCLINIC INTEGRATION
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use dyn_mod
use work_mod
use pmix_mod
use msg_mod
!lhl1204
use forc_mod, only: psa,USTAR,BUOYTUR, BUOYSOL,NSWV,SWV
!lhl1204
 
      IMPLICIT NONE
      REAL(r8)   :: ABCD,TUP,SUP,TLO,SLO,RHOUP,RHOLO
      REAL(r8)   :: DENS
      real(r8),dimension(imt,jmt,km)::ALPHA,BETA,pp
      real(r8),dimension(imt,jmt,km)::ppa,ppb,ppc
      EXTERNAL DENS
 
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            H0L (I,J)= H0F (I,J)
            H0F (I,J)= H0 (I,J)
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UTL (I,J,K)= UTF (I,J,K)
               VTL (I,J,K)= VTF (I,J,K)
               UTF (I,J,K)= U (I,J,K)
               VTF (I,J,K)= V (I,J,K)
            END DO
         END DO
      END DO
 
!     --------------------------------------------------------------
!     PREPARING FOR CALCULATING RICHARDSON NUMBER
!     --------------------------------------------------------------
 
      allocate(dlu(imt,jmt,km),dlv(imt,jmt,km),gg(imt,jmt,km))
      allocate(rit(imt,jmt,kmm1),ric(imt,jmt,kmm1))
      allocate(rict(imt,jmt,kmm1))
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JEM
            DO I = 2,IMM
               DLU (I,J,K)= 0.25D0* (AT (I,J,K,1) + AT (I -1,J,K,1) + &
                                AT (I,J +1,K,1) +  AT (I -1,J +1,K,1))*VIV(I,J,K)
               DLV (I,J,k)= 0.25D0* (AT (I,J,K,2) + AT (I -1,J,K,2) + &
                               AT ( I,J +1,K,2) +  AT (I -1,J +1,K,2))*VIV(I,J,K)
!lhl0608               DLU (I,J,K)= 0.25D0* (ATB(I,J,K,1) + ATB(I -1,J,K,1) + &
!lhl0608                                ATB(I,J +1,K,1) +  ATB(I -1,J +1,K,1))*VIV(I,J,K)
!lhl0608               DLV (I,J,k)= 0.25D0* (ATB(I,J,K,2) + ATB(I -1,J,K,2) + &
!lhl0608                               ATB( I,J +1,K,2) +  ATB(I -1,J +1,K,2))*VIV(I,J,K)
            END DO
!               if (nx_proc == 1) then
!               DLU (1,J,K)= DLU (IMM,J,K)
!               DLU (IMT,J,K)= DLU (2,J,K)
!               DLV (1,J,K)= DLV (IMM,J,K)
!               DLV (IMT,J,K)= DLV (2,J,K)
!               end if
         END DO
!lhl0710            DO I = 1,IMT
!lhl0710               DLU (I,JMT,K)= DLU (I,JEM,K)
!lhl0710               DLV (I,JMT,K)= DLV (I,JEM,K)
!lhl0710            END DO
      END DO

     if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (K,J)
           do k=1,km
           do j=jst,jem
               DLU (1,J,K)= DLU (IMM,J,K)
               DLU (IMT,J,K)= DLU (2,J,K)
               DLV (1,J,K)= DLV (IMM,J,K)
               DLV (IMT,J,K)= DLV (2,J,K)
           end do
           end do
      end if

#ifdef SPMD
!     call mpi_barrier(mpi_comm_ocn,ierr)
      call exch_boundary(DLU,km)
      call exch_boundary(DLV,km)
#endif
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KMM1
         DO J = JST,JET
            DO I = 1,IMT
               rit (i,j,k)= 0.0D0
               ric (i,j,k)= 0.0D0 ! Dec. 18, 2002, Yongqiang YU!
!lhl1204
               rict (i,j,k)= 0.0D0
               ricdt (i,j,k)= 0.0D0
!lhl1204
            END DO
         END DO
      END DO

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               akt (i,j,k,1)= 0.0D0
               akt (i,j,k,2)= 0.0D0
            END DO
         END DO
      END DO
 
!
! ric locate at integer level
!

      DO K = 1,KMM1
!         DO J = JSM,JEM ! Dec. 11, 2002, Yongqiang YU
!$OMP PARALLEL DO PRIVATE(J,I,tup,sup,TLO,SLO,rhoup,RHOLO), shared(k,dlu,dlv,to,so,ric,od0,g,odzt)
         DO J = JST,JET
            DO I = 2,IMM
               TUP = DLU (I,J,K) - TO (K +1)
               SUP = DLV (I,J,K) - SO (K +1)
               TLO = DLU (I,J,K +1) - TO (k +1)
               SLO = DLV (I,J,K +1) - SO (K +1)
               RHOUP = DENS (TUP, SUP, K +1)
               RHOLO = DENS (TLO, SLO, K +1)
!lhl1204
!               ric (I,J,K) = VIV (I,J,K +1)* OD0* G * (RHOLO - RHOUP)
               ric (I,J,K) = VIV (I,J,K +1)* OD0* G * (RHOLO - RHOUP)*ODZT(K+1)
!lhl1204
            END DO
         END DO
      END DO
!
      if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE(j,k)
        do k=1,kmm1
        do j=jst,jmt
          ric(1,j,k)=ric(imm,j,k)
          ric(imt,j,k)=ric(2,j,k)
        end do
        end do
      end if
! 
#ifdef SPMD
      call exch_boundary(ric,kmm1)
#endif
#ifdef DEBUG
      call chk_var1(ric,'ic',0)
#endif

!lhl1204
!
! ric locate at integer level from 2
!
      DO K = 1,KMM1
!$OMP PARALLEL DO PRIVATE (J,I,TUP,SUP,TLO,SLO,RHOUP,RHOLO), shared(k,to,so,at,rict,OD0,G)
         DO J = JST,JET
            DO I = 2,IMM
               TUP = AT (I,J,K,1) - TO (K +1)
               SUP = AT (I,J,K,2) - SO (K +1)
               TLO = AT (I,J,K +1,1) - TO (k +1)
               SLO = AT (I,J,K +1,2) - SO (K +1)
!lhl0608               TUP = ATB(I,J,K,1) - TO (K +1)
!lhl0608               SUP = ATB(I,J,K,2) - SO (K +1)
!lhl0608               TLO = ATB(I,J,K +1,1) - TO (k +1)
!lhl0608               SLO = ATB(I,J,K +1,2) - SO (K +1)
               RHOUP = DENS (TUP, SUP, K +1)
               RHOLO = DENS (TLO, SLO, K +1)
               rict (I,J,K) = VIT (I,J,K +1)* OD0* G * (RHOLO - RHOUP)*ODZT(K+1)
            END DO
         END DO
      END DO
!
     if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE(j,k)
     do k=1,kmm1
     do j = jst, jmt
         rict(1,j,k)=rict(imm,j,k)
         rict(imt,j,k)=rict(2,j,k)
     end do
     end do
     end if
#ifdef SPMD
      call exch_boundary(rict,kmm1)
#endif
#ifdef DEBUG
      call chk_var1(rict,'ct',1)
#endif
!lhl1204
 
!     --------------------------------------------------------------
!     COMPUTING DENSITY AND BAROCLINIC PRESSURE
!     --------------------------------------------------------------
 
! calculate the potential density using BC equation on U-grid
! note the referrence of density was added!
! if the salinity should multiple 1000??????
      CALL DENSITY

!lhl0711#ifdef SPMD
!lhl0711      call exch_boundary(PDENSITY,km)
!lhl0711#endif

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
!lhl1204
!               GG (I,J,K)= - OD0* G * WKA (I,J,K)
               GG (I,J,K)=-OD0*G*PDENSITY (I,J,K)*VIT(I,J,K)
!lhl1204
            END DO
         END DO
      END DO

!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
!lhl1204
! calculate pressure at level 1
!            WKA (I,J,1)= GG (I,J,1)*0.5D0* DZP (1)* VIT (I,J,1)
            PP (I,J,1)= GG (I,J,1)*0.5D0* DZP (1)* VIT (I,J,1)
!
            PPA(I,J,1)= PSA (I,J)* VIT (I,J,1)
            PPB(I,J,1)= AT(I,J,1,1)*VIT(I,J,1)
            PPC(I,J,1)= AT(I,J,1,2)*VIT(I,J,1)
!lhl0608            PPB(I,J,1)= ATB(I,J,1,1)*VIT(I,J,1)
!lhl0608            PPC(I,J,1)= ATB(I,J,1,2)*VIT(I,J,1)
!lhl1204
         END DO
      END DO


      DO K = 2,KM
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JST,JET
            DO I = 1,IMT
!lhl1204
! calculate pressure at T-grid
!               WKA (I,J,K)= VIT (I,J,K)* (WKA (I,J,K -1) +0.5D0* &
!               (GG (I,J,K)* DZP (K) + GG (I,J,K -1)* DZP (K -1)))
               PP (I,J,K)= VIT (I,J,K)* (PP(I,J,K -1) +0.5D0* &
               (GG (I,J,K)* DZP (K) + GG (I,J,K -1)* DZP (K -1)))
!
               PPA(I,J,K)= VIT (I,J,K)* (PPA(I,J,K -1)+GG (I,J,K -1)* DZP (K -1))
               PPB(I,J,K)= VIT (I,J,K)* (AT (I,J,K-1,1)+(AT (I,J,K-1,1)-AT (I,J,K,1))* &
                               DZP (K -1)/(DZP(K-1)+DZP(K)))
               PPC(I,J,K)= VIT (I,J,K)* (AT (I,J,K-1,2)+(AT (I,J,K-1,2)-AT (I,J,K,2))* &
                               DZP (K -1)/(DZP(K-1)+DZP(K)))
!lhl0608               PPB(I,J,K)= VIT (I,J,K)* (ATB(I,J,K-1,1)+(ATB(I,J,K-1,1)-ATB(I,J,K,1))* &
!lhl0608                               DZP (K -1)/(DZP(K-1)+DZP(K)))
!lhl0608               PPC(I,J,K)= VIT (I,J,K)* (ATB(I,J,K-1,2)+(ATB(I,J,K-1,2)-ATB(I,J,K,2))* &
!lhl0608                               DZP (K -1)/(DZP(K-1)+DZP(K)))
!lhl1204
            END DO
         END DO
      END DO
!lhl0711#ifdef SPMD
!lhl0711      call exch_boundary(PP,km)
!lhl0711      call exch_boundary(PPA,km)
!lhl0711      call exch_boundary(PPB,km)
!lhl0711      call exch_boundary(PPC,km)
!lhl0711#endif
#ifdef DEBUG
      call chk_var3d(pp,'pp',1)
      call chk_var3d(ppa,'pa',1)
      call chk_var3d(ppb,'pb',1)
      call chk_var3d(ppc,'pc',1)
#endif
!lhl1204
! calculate the thermal expansion and the salinity contraction at T-grid
!  PP in negative in model
!  in PP a OD0 is mutilped and should be divided
!  PP is in Pa to tranform to db by divided 10000.
!
!      CALL THERMAL(AT(1,1,1,1),AT(1,1,1,2),PP,ALPHA,BETA,VIT)
      CALL THERMAL(PPB,PPC,PPA,ALPHA,BETA,VIT)
!
!lhl0711#ifdef SPMD
!lhl0711      call exch_boundary(ALPHA,km)
!lhl0711      call exch_boundary(BETA,km)
!lhl0711#endif
#ifdef DEBUG
      call chk_var3d(alpha,"al",1)
      call chk_var3d(beta,"be",1)
#endif
!
! calculate the surface buoyancy fluxes at T-grid

!M
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
         BUOYTUR(I,J)= VIT(I,J,1)*NSWV(I,J)*G*ALPHA(I,J,1)*OD0CP
         BUOYSOL(I,J)= VIT(I,J,1)* SWV(I,J)*G*ALPHA(I,J,1)*OD0CP
         END DO
      END DO
!
#ifdef DEBUG
      call chk_var2d(BUOYTUR,"ur",1)
      call chk_var2d(BUOYSOL,"ol",1)
#endif
!
! calculate the T component minus S component of B-V frequency at T-grid
!

!M
!$OMP PARALLEL DO PRIVATE (k,J,I)
      DO K = 1,KMM1
         DO J = JST,JET
            DO I = 1,IMT
               ricdt(I,J,K)=VIT(I,J,K+1)*G*((AT(I,J,K,1)-AT(I,J,K+1,1))*ALPHA(I,J,K+1)&
                                     +1000.D0*(AT(I,J,K,2)-AT(I,J,K+1,2))*BETA(I,J,K+1))*ODZT(K+1)
!lhl0608               ricdt(I,J,K)=VIT(I,J,K+1)*G*((ATB(I,J,K,1)-ATB(I,J,K+1,1))*ALPHA(I,J,K+1)&
!lhl0608                                     +1000.D0*(ATB(I,J,K,2)-ATB(I,J,K+1,2))*BETA(I,J,K+1))*ODZT(K+1)
!
!               rict(I,J,K)=VIT(I,J,K+1)*G*((AT(I,J,K,1)-AT(I,J,K+1,1))*ALPHA(I,J,K+1)&
!                                     -1000.D0*(AT(I,J,K,2)-AT(I,J,K+1,2))*BETA(I,J,K+1))*ODZT(K+1)
            END DO
         END DO
      END DO
!lhl0711#ifdef SPMD
!lhl0711      call exch_boundary(ricdt,kmm1)
!lhl0711!      call exch_boundary(rict,kmm1)
!lhl0711#endif
#ifdef DEBUG
      call chk_var1(ricdt,"dt",1)
!      call chk_var1(rict,"ct",1)
#endif

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               GG (I,J,K)=-OD0*G*(PDENSITY (I,J,K)-PO(K)-1000.0D0)*VIT(I,J,K)
            END DO
         END DO
      END DO
!lhl1204
 
!     --------------------------------------------------------------
!     COMPUTING HORIZONTAL GRADIENT OF BAROCLINIC PRESSURE
!     --------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMT
!lhl1204
!               DLV (I,J,K)= OUY (J)*0.5D0* &
!               (WKA (I,J +1,K) - WKA (I,J,K) + WKA (I -1,J +1,K) - WKA (I -1,J,K))
!               DLU (I,J,K)= OUX (J)*0.5D0* &
!               (WKA (I,J,K) - WKA (I -1,J,K) + WKA (I,J +1,K) - WKA (I -1,J +1,K))
               DLV (I,J,K)= OUY (J)*0.5D0* &
               (PP (I,J +1,K) - PP (I,J,K) + PP (I -1,J +1,K) - PP (I -1,J,K))
               DLU (I,J,K)= OUX (J)*0.5D0* &
               (PP (I,J,K) - PP (I -1,J,K) + PP (I,J +1,K) - PP (I -1,J +1,K))
!lhl1204
            END DO
         END DO
      END DO
 
     if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (K,J)
        do J=JSM,JEM
        DO K=1,KM
            DLU (1,J,K)= DLU (IMM,J,K)
            DLV (1,J,K)= DLV (IMM,J,K)
        end do
        end do
     end if
#ifdef  SPMD                                                   
       call exch_boundary(dlu,km)
       call exch_boundary(dlv,km)
#endif 
 
!     --------------------------------------------------------------
!     COMPUTING VERTICALLY INTEGRATED GRADIENT OF BAROCLINIC PRESSURE
!     --------------------------------------------------------------
 
      CALL VINTEG (DLU,PXB)
      CALL VINTEG (DLV,PYB)

 
!     --------------------------------------------------------------
!     COMPUTING WGP, WHX, WHY, PAX & PAY
!     --------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            DLU (I,J,1)= 0.0D0
            DLU (I,J,2)= 0.0D0
         END DO
      END DO
 
      DO K = 1,KM
!$OMP PARALLEL DO PRIVATE (J,I,ABCD)
        DO J = JST,JET
         DO I = 1,IMT
               ABCD = GG (I,J,K)* OHBT (I,J)* DZP (K)
               DLU (I,J,1)= DLU (I,J,1) + ABCD
               DLU (I,J,2)= DLU (I,J,2) + ABCD * ZKT (K)
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            DLV (I,J,1)= (DLU (I,J,1) + DLU (I,J,2)* OHBT (I,J))/ G
            DLV (I,J,2)= DLU (I,J,2)* OHBT (I,J)* OHBT (I,J)
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JSM,JEM
         DO I = 2,IMM
            WGP (I,J)= 0.25D0* VIV (I,J,1)* &
            (DLV (I,J,1) + DLV (I -1,J,1) + DLV (I,J +1,1) + DLV (I -1, &
             J +1,1))     
         END DO
      END DO

     if (nx_proc == 1) then
     do j=jsm,jem
         WGP (1,J)= WGP (IMM,J)
         WGP (IMT,J)= WGP (2,J)
     end do
     end if
 
!$OMP PARALLEL DO PRIVATE (J,I,ABCD)
      DO J = JSM,JEM
         DO I = 2,IMM
            ABCD = 0.25D0* VIV (I,J,1)* &
                   (DLV (I,J,2) + DLV (I -1,J,2) + DLV (I,J +1,2)       &
                     + DLV (I -1,J +1,2))
            WHX (I,J)= HBX (I,J)* ABCD
            WHY (I,J)= HBY (I,J)* ABCD
         END DO
      END DO
!
    if (nx_proc == 1) then
    do j=jsm,jem
         WHX (1,J)= WHX (IMM,J)
         WHX (IMT,J)= WHX (2,J)
         WHY (1,J)= WHY (IMM,J)
         WHY (IMT,J)= WHY (2,J)
    end do
    end if

#ifdef  SPMD                                                   
       call exch_boundary(WGP,1)
       call exch_boundary(WHX,1)
       call exch_boundary(WHY,1)
#endif  
 
!$OMP PARALLEL DO PRIVATE (J,I),shared(OD0)
      DO J = JSM,JEM
         DO I = 2,IMM
            PAY (I,J)= - OD0* OUY (J)*0.5D0* &
            (PSA (I,J +1) - PSA (I,J) + PSA (I -1,J +1) - PSA (I -1,J)) 
            PAX (I,J)= - OD0* OUX (J)*0.5D0* &
            (PSA (I,J) - PSA (I -1,J) + PSA (I,J +1) - PSA (I -1,J +1))
         END DO
      END DO

      RETURN
      END SUBROUTINE READYT
 
 
