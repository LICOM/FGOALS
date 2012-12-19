!     =================
      SUBROUTINE READYC
!     =================
!     PREPARATION OF BAROTROPIC AND BAROCLINIC INTEGRATION
 
!     ADVECTION + DIFFUSION
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use tracer_mod
use pmix_mod
use forc_mod, only: su,sv,USTAR,BUOYTUR, BUOYSOL
 
      IMPLICIT NONE
!      REAL(r8)  :: WKP (KMP1)
      INTEGER   :: IWK,n2
      REAL(r8)  :: WK1 (KM) ,WK2 (KM), WK3 (KM)
      REAL(r8)  :: WP1 (KM) ,WP2 (KM), WP3 (KM)
      REAL(r8)  :: WP4 (KM) ,WP5 (KM), WP6 (KM)
      REAL(r8)  :: WP7 (KM) ,WP8(KM*2)
      REAL(r8)  :: WP9 ,WP10, WP11
      REAL(r8),dimension(IMT,JMT,KM) :: WP12,WP13
      REAL(r8)  :: riv1,riv2,epsln,RKV,RKV1
      REAL(r8)  :: adv_x1,adv_x2,adv_x,adv_y1,adv_y2,adv_z,diff_u1,diff_u2,diff_v1,diff_v2
      REAL(r8)  :: dlux,dlvx,dluy,dlvy,dluz,dlvz,adv_z1,adv_z2,adv_z3,adv_z4
      REAL(r6)  :: xxx
      
       
#if (defined CANUTO)
      REAL(r8)  :: AIDIF
#endif

!YU 
      real (r8):: akt_back(2*km),aks_back(2*km),akm_back(2*km)
      allocate(riu(imt,jmt,0:km),stat=ierr)
!Yu   if (ierr /= 0) then
!Yu      write(6,*)'allocation error---riu'
!Yu      stop
!Yu   end if

 
#if (defined BIHAR)
      allocate(tmp1(imt,jmt,km),tmp2(imt,jmt,km))
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               TMP1 (I,J,K)= 0.0D0
               TMP2 (I,J,K)= 0.0D0
            END DO
         END DO
      END DO
 
#endif
 
      epsln = 1.0D-25
      
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            H0BL (I,J)= H0BF (I,J)
            H0BF (I,J)= H0 (I,J)
         END DO
      END DO
 
!lhl0711
#if (defined CANUTO)
       AIDIF=0.0  
!      if (ISC/=0)  AIDIF=0.5  
#endif
!lhl0711
 
!---------------------------------------------------------------------
!     Calculating Richardson number riu (at U/V-point);
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 1,kmm1
         DO J = jst,jet 
            DO I = 1,imt
            s2t(i,j,k) = 0.D0
            ridt(i,j,k) = 0.D0
            END DO
         END DO
         END DO

!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 0,km
         DO J = jst,jet
            DO I = 1,imt
            riu (i,j,k) = 0.D0
            END DO
         END DO   
         END DO

!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 1,km
         DO J = jst,jet
            DO I = 1,imt
            wp12(i,j,k) = 0.D0
            wp13(i,j,k) = 0.D0
            END DO
         END DO   
         END DO

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               wp12(i,j,k)=   VIT(I  ,J,K)*&
               (up(i,j,k  )+up(i+1,j,k  )+up(i,j-1,k  )+up(i+1,j-1,k  )) &
              /(VIV(i,j,k)+VIV(i+1,j,k)+VIV(i,j-1,k)+VIV(i+1,j-1,k)+epsln)
               wp13(i,j,k)=   VIT(I  ,J,K)*&
               (vp(i,j,k  )+vp(i+1,j,k  )+vp(i,j-1,k  )+vp(i+1,j-1,k  )) &
              /(VIV(i,j,k)+VIV(i+1,j,k)+VIV(i,j-1,k)+VIV(i+1,j-1,k)+epsln)
            END DO
         END DO
      END DO
!
      if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (K,J)
         do k=1,km
         do j=jsm,jem
              wp12(1,j,k)=wp12(imm,j,k)
              wp12(imt,j,k)=wp12(2,j,k)
              wp13(1,j,k)=wp13(imm,j,k)
              wp13(imt,j,k)=wp13(2,j,k)
         end do
         end do
      end if
!lhl0711#ifdef SPMD 
!lhl0711      call exch_boundary(wp12,km)
!lhl0711      call exch_boundary(wp13,km)
!lhl0711#endif      
 
!$OMP PARALLEL DO PRIVATE (K,J,I,riv1,riv2)
      DO K = 1,KMM1
         DO J = JST,JET
            DO I = 2,IMM
               riv1 = wp12 (I,J,K) - wp12 (I,J,K +1)
               riv2 = wp13 (I,J,K) - wp13 (I,J,K +1)
               s2t (i,j,k) =vit(i,j,k+1)*(riv1*riv1+riv2*riv2)*ODZT(K+1)*ODZT(K+1)
               ridt(i,j,k) =vit(i,j,k+1)*ricdt(i,j,k)/(s2t(i,j,k)+epsln)
#ifdef CANUTO  
               rit (i,j,k)= VIT (I,J,K +1)*rict(i,j,k)/(s2t(i,j,k)+epsln)
#else          
               rit (i,j,k)= rit (i,j,k) +VIT (I,J,K +1)*rict(i,j,k)/(s2t(i,j,k)+epsln)
#endif
            END DO
         END DO
      END DO
      if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (K,J)
         do k=1,kmm1
         do j=jst,jet
            s2t(1,j,k) = s2t(IMM,j,k)
            s2t(imt,j,k) = s2t(2,j,k)
            rit(1,j,k)=rit(imm,j,k)
            rit(imt,j,k)=rit(2,j,k)
            ridt(1,j,k)=ridt(imm,j,k)
            ridt(imt,j,k)=ridt(2,j,k)
         end do
         end do
      endif
!
!lhl0711#ifdef SPMD
!lhl0711      call exch_boundary(s2t,kmm1)
!lhl0711      call exch_boundary(rit,kmm1)
!lhl0711      call exch_boundary(ridt,kmm1)
!lhl0711#endif
#ifdef DEBUG
      call chk_var1(s2t,'2t',1)
      call chk_var1(rit,'it',1)
      call chk_var1(ridt,'dt',1)
#endif
!
!$OMP PARALLEL DO PRIVATE (K,J,I,riv1,riv2)
      DO K = 1,KMM1
!Yu      DO J = JSM,JEM
         DO J = JST,JET    ! Dec. 4, 2002, Yongqiang Yu
            DO I = 2,IMM
               riv1 = UP (I,J,K) - UP (I,J,K +1)
               riv2 = VP (I,J,K) - VP (I,J,K +1)
               s2u (i,j,k) =viv(i,j,k+1)*(riv1*riv1+riv2*riv2)*ODZT(K+1)*ODZT(K+1)
! calculate the shear square and the T component minus S component of Richardson Number
!lhl1204
               riu (i,j,k) = VIV (I,J,K +1)*ric (i,j,k)/(s2u(i,j,k)+epsln)
!               riu (i,j,k) = ric (i,j,k)/(riv1*riv1/ODZT(K+1)+riv2*riv2/ODZT(K+1)+epsln)
!lhl1204
            END DO
            riu (1,j,k) = riu (IMM,j,k)
            riu (imt,j,k) = riu (2,j,k)
         END DO
      END DO
!
!lhl0711#ifdef SPMD
!lhl0711      call exch_boundary(riu,kmm1)
!lhl0711#endif
#ifdef DEBUG
      call chk_var1(riu,'iu',0)
#endif
!
!         DO J = 1,jmt
!         DO I = 2,IMM
!      if (VIT(I,J,1).gt.0.5) write(*,'(2i4,11f10.3)')i,j,vit(i,j,1),(rit(i,j,k),k=1,10)
!      if (VIT(I,J,1).gt.0.5) write(*,'(2i4,11f10.3)')i,j,vit(i,j,1),(ridt(i,j,k),k=1,10)
!      if (VIT(I,J,1).gt.0.5) write(*,'(2i4,11f10.3)')i,j,vit(i,j,1),(s2t(i,j,k)*1d+6,k=1,10)
!      if (VIT(I,J,1).gt.0.5) write(*,'(2i4,11f10.3)')i,j,vit(i,j,1),(rict(i,j,k)*1d+6,k=1,10)
!      if (VIT(I,J,1).gt.0.5) write(*,'(2i4,11f10.3)')i,j,vit(i,j,1),(rict(i,j,k)/s2t(i,j,k),k=1,10)
!         END DO
!         END DO
!
#ifdef CANUTO
      DO J = JST,JET
         DO I = 1,IMT
         AMLD(I,J)=0.0D0
         DO K = 1,KM
         AKMT(I,J,K)=0.0D0
         AKMU(I,J,K)=0.0D0
!         AKT(I,J,K,1)=0.0D0
!         AKT(I,J,K,2)=0.0D0
         END DO
         END DO
      END DO
!
!
#ifdef DEBUG
      call chk_var3d(at(1,1,1,1),"tt",1)
      call chk_var3d(at(1,1,1,2),"ss",1)
      call chk_var3d(pdensity,"pd",1)
      call chk_var3d(up,"us",0)
      call chk_var3d(vp,"vs",0)
!
      call chk_var1(rit,"it",1)
      call chk_var1(ridt,"dt",1)
      call chk_var1(s2t,"2t",1)
      call chk_var1(rict,"ct",1)
!
      call chk_var2d(ustar,"us",1)
      call chk_var2d(buoytur,"br",1)
      call chk_var2d(buoysol,"bs",1)
#endif
!
!         DO K = 1,ITNU(55,4)-1
!         print*,k
!         print*,vit(55,4,k)
!         print*,RIT(55,4,k)
!         print*,RIDT(55,4,k)
!         print*,S2T(55,4,k)
!         print*,RICT(55,4,k)
!         print*
!         END DO
!lhl241204

!$OMP PARALLEL DO PRIVATE (J,I,K,wp1,wp2,wp3,wp4,wp5,wp6,wp7,wp8,wp9,wp10,wp11,akt_back,akm_back,aks_back,iwk,wk1,wk2,wk3,xxx)
      DO J = JSM,JEM
         DO I = 2,IMM
         !      DO J = 80,90
!         DO I = 180,200
!        if (ITNU(I,J).gt.0.8) then
        if (VIT(I,J,1).gt.0.5) then
!
         wp1=0.D0
         wp2=0.D0
         wp3=0.D0
         wp4=0.D0
         wp5=0.D0
         wp6=0.D0
         wp7=0.D0
         wp8=0.D0
         wp9=0.D0
         wp10=0.D0
         wp11=0.D0
!
         do k=1,km
            AKM_BACK(K)=0.0D0
            AKT_BACK(K)=0.0D0
            AKS_BACK(K)=0.0D0
         end do
!
         DO K = 1,ITNU(i,j)-1
               wp1(K)= VIT (I,J,K+1)* (AT (I,J,K,1)+(AT (I,J,K,1)-AT (I,J,K+1,1))* &
                               DZP (K)/(DZP(K)+DZP(K+1)))
               wp2(K)= VIT (I,J,K+1)* (AT (I,J,K,2)+(AT (I,J,K,2)-AT (I,J,K+1,2))* &
                               DZP (K)/(DZP(K)+DZP(K+1)))
!lhl0608               wp1(K)= VIT (I,J,K+1)* (ATB(I,J,K,1)+(ATB(I,J,K,1)-ATB(I,J,K+1,1))* &
!lhl0608                               DZP (K)/(DZP(K)+DZP(K+1)))
!lhl0608               wp2(K)= VIT (I,J,K+1)* (ATB(I,J,K,2)+(ATB(I,J,K,2)-ATB(I,J,K+1,2))* &
!lhl0608                               DZP (K)/(DZP(K)+DZP(K+1)))
               wp3(K)= VIT (I,J,K+1)* (pdensity (I,J,K)+(pdensity (I,J,K)-pdensity(I,J,K+1))* &
                               DZP (K)/(DZP(K)+DZP(K+1)))*1.d-3
               wp8(k)=-vit(i,j,k+1)*ZKP(k+1)*1.d+2
!         wp1(k)=vit(i,j,k)*at(i,j,k,1)
!         wp2(k)=vit(i,j,k)*at(i,j,k,2)
!         wp3(k)=vit(i,j,k)*pdensity(i,j,k)*1.d-3
!         wp8(k)=-vit(i,j,k)*ZKT(k)*1.d+2
         END DO
         DO K = 1,ITNU(i,j)-1
         wp4(k)=vit(i,j,k+1)*RIT(i,j,k)
         wp5(k)=vit(i,j,k+1)*RIDT(i,j,k)
         wp6(k)=vit(i,j,k+1)*S2T(i,j,k)
         wp7(k)=vit(i,j,k+1)*RICT(i,j,k)
         END DO
         wp9=vit(i,j,1)*USTAR(I,J)*1.0d+2
         wp10=vit(i,j,1)*BUOYTUR(I,J)*1.0d+4
         wp11=vit(i,j,1)*BUOYSOL(I,J)*1.0d+4
!
         IWK=ITNU(I,J)-1
!
!     if (mytid.eq.0) then
!       print*,'out'
!       print*,i,j 
!       print*,dzp
!     endif
!input ZKT in cm, AT(1) in C, AT(2) in (s-35)/1000., PDENSITY in g/cm^3


         CALL  TURB_2(wp8,wp1,wp2,wp3,&
!input RIT and RIDT no unit, S2T in 1/s^2, DFRICMX and DWNDMIX in cm^2/s
                wp4,wp5,wp6, DFRICMX*1.0d+4,DWNDMIX*1.0d+4,&
!output in cm^2/s, so 1d-4 should be multipled
               AKM_BACK,AKT_BACK,AKS_BACK,&
!input  RICT in 1/s^2 USTAR in cm/s, BUOYTUR,BUOYSOL in cm^2/s^3,FF in 1/s
               wp7,wp9,wp10,wp11,FF(J),& !OK
!output amld in cm, akmt,akh, and aks in cm^2/s
!               AMLD(I,J),AKMT(I,J,1),AKT(I,J,1,1),AKT(I,J,1,2),&
               AMLD(I,J),WK1,WK2,WK3,&
!input int
               IWK,NA(I,J),KM,1,0,0,i,j) !OK!
             
!           if ( i==8 .and. j==4) then
!           write(*,*) WK1
!           write(*,*) WK2
!           write(*,*) WK3
!           stop
!           end if
      


         DO K = 1,KM
!
         xxx = WK1(K)
         WK1(K) = xxx
         xxx = WK2(K)
         WK2(K) = xxx
         xxx = WK3(K)
         WK3(K) = xxx
!
!        if ( k < 10) then
!        if ( j_global(j) > 1160 .and. j_global(j) < 1280 .and. i_global(i) > 1480 .and. i_global(i) < 1560 ) then
!           wk1(k)= max(wk1(k),0.1D+2)
!           wk2(k)= max(wk2(k),0.1D+2)
!           wk2(k)= max(wk2(k),0.1D+2)
!        end if
!        end if
!         AKMT(I,J,K)=+WK1(K)*1.d-4
!         AKT(I,J,K,1)=AKT(I,J,K,1)+WK2(K)/NCC*1.d-4
!         AKT(I,J,K,2)=AKT(I,J,K,2)+WK3(K)/NCC*1.d-4
         AKMT(I,J,K)=+(WK1(K)+dmin1(AKM_BACK(K),1d-3))*1.d-4
         AKT(I,J,K,1)=AKT(I,J,K,1)+(WK2(K)+dmin1(AKT_BACK(K),1d-3))/NCC*1.d-4
         AKT(I,J,K,2)=AKT(I,J,K,2)+(WK3(K)+dmin1(AKS_BACK(K),1d-3))/NCC*1.d-4

!         AKMT(I,J,K)=+WK1(K)+AMV
!         AKT(I,J,K,1)=AKT(I,J,K,1)+WK2(K)+AHV
!         AKT(I,J,K,2)=AKT(I,J,K,2)+WK3(K)+AHV
         END DO
!       else
!         AMLD(I,J)=0.0D0
!         DO K = 1,KM
!         AKM_BACK(I,J,K)=0.0D0
!         AKT_BACK(I,J,K)=0.0D0
!         AKS_BACK(I,J,K)=0.0D0
!         AKMT(I,J,K)=0.0D0
!         AKMU(I,J,K)=0.0D0
!         AKT(I,J,K)=0.0D0
!         AKS(I,J,K)=0.0D0
!         END DO
!     if (mytid.eq.0) then
!       print*,'out'
!       print*,i,j
!       write(*,'(10d15.5)')(wk2(k),k=1,na(i,j))
!     endif
        endif
         END DO
      END DO

     if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (J)
         do j=jsm,jem
              amld(1,j)=amld(imm,j)
              amld(imt,j)=amld(2,j)
         end do
!$OMP PARALLEL DO PRIVATE (K,J)
         DO K = 1,KM
         do j=jsm,jem
              akmt(1,j,k)=akmt(imm,j,k)
              akmt(imt,j,k)=akmt(2,j,k)
              akt(1,j,k,1)=akt(imm,j,k,1)
              akt(imt,j,k,1)=akt(2,j,k,1)
              akt(1,j,k,2)=akt(imm,j,k,2)
              akt(imt,j,k,2)=akt(2,j,k,2)
         END DO
         END DO
     end if

#ifdef SPMD
      call exch_boundary(amld,1)
      call exch_boundary(akmt,km)
      call exch_boundary(akt(1,1,1,1),km)
      call exch_boundary(akt(1,1,1,2),km)
#endif


!!
!     do j=jsm,jem
!     write(125,*) "OK---",j, (akmt(i,j,5),i=1,200)
!     end do
!     close(125)
!     stop
!!
! calculate the vertical mixing on U-grid
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KMM1
         DO J = JSM,JEM
            DO I = 2,IMM
               AKMU (I,J,K)= 0.25D0* (AKMT (I,J,K) + AKMT (I -1,J,K)&
                    + AKMT (I,J +1,K) + AKMT (I -1,J +1,K))*VIV(I,J,K+1)
            END DO
         END DO
      END DO
      if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (K,J)
           do k=1,kmm1
           do j=jsm,jem
              AKMU(1,j,k)=AKMU(IMM,j,k)
              AKMU(IMT,j,k)=AKMU(2,j,k)
           end do
           end do
      end if
#ifdef SPMD
      call exch_boundary(akmu(1,1,1),km)
#endif
#ifdef DEBUG
      call chk_var3d(akmt,"km",1)
      call chk_var3d(akmu,"km",0)
      call chk_var3d(akt(1,1,1,1),"kt",1)
      call chk_var3d(akt(1,1,1,2),"ks",1)
      call chk_var2d(amld,"md",1)
#endif
!lhl241204
#endif
!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM: ZONAL COMPONENT
!---------------------------------------------------------------------
      CALL UPWELL (U,V,H0)
      
 
 
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
 
!Yu   if (ierr /= 0) then
!Yu      write(6,*)'allocation error---dlu,dlv'
!Yu      stop
!Yu   end if

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JET
            DO I = 1,IMT
               DLU (I,J,K)= 0.0D0
               DLV (I,J,K)= 0.0D0
               WKA (I,J,K)= 0.0D0
            END DO
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,K)= 0.25D0* (WS (I,J,K) + WS(I-1,J,K) + WS(I,J+1,K)+&
               WS (I -1,J +1,K))             
            END DO
         END DO
      END DO
 

!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM: ZONAL COMPONENT
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I,adv_x1,adv_x2,adv_y1,adv_y2,adv_z1,adv_z2,adv_z3,adv_z4, &
!$OMP                      dlux,dlvx,dluy,dlvy, dluz,dlvz)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               adv_x1=0.25D0* OUX (J)* (U (I -1,J,K) + U (I,J,K))
               adv_x2=0.25D0* OUX (J)* (U (I +1,J,K) + U (I,J,K))
               dlux = - ( adv_x2* (U (I +1,J,K) - U (I,J,K)) + &
                                 adv_x1* (U (I,J,K) - U (I -1,J,K))) 
               dlvx = - ( adv_x2* (V (I +1,J,K) - V (I,J,K)) + &
                                 adv_x1* (V (I,J,K) - V (I -1,J,K)))
               adv_y1=V (I,J -1,K) + V (I,J,K)
               adv_y2=V (I,J +1,K) + V (I,J,K)
               dluy = -( R1B (J)* adv_y2* (U (I,J +1,K) - U (I,J,K)) + &
                             R1A (J)* adv_y1 * (U (I,J,K) - U (I,J -1,K)))
               dlvy = -( R1B (J)* adv_y2* (V (I,J +1,K) - V (I,J,K)) + &
                             R1A (J)* adv_y1 * (V (I,J,K) - V (I,J -1,K)))
               if (k==1 )then
                   adv_z1=0.0D0
                   adv_z3=0.0D0
               else
                   adv_z1=WKA (I,J,K)* (U (I,J,K -1) - U (I,J,K))
                   adv_z3=WKA (I,J,K)* (V (I,J,K -1) - V (I,J,K))
               end if
               if (k==km )then
                   adv_z2=0.0D0
                   adv_z4=0.0D0
               else
                   adv_z2=WKA (I,J,K+1)* (U (I,J,K ) - U (I,J,K+1))
                   adv_z4=WKA (I,J,K+1)* (V (I,J,K ) - V (I,J,K+1))
               end if
               dluz = - 0.5D0* ODZP (K)* (adv_z1+adv_z2)
               dlvz = - 0.5D0* ODZP (K)* (adv_z3+adv_z4)
               dlu(i,j,k)=dlux+dluy+dluz
               dlv(i,j,k)=dlvx+dlvy+dlvz
            END DO
         END DO
      END DO
!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM: ZONAL COMPONENT
!---------------------------------------------------------------------
 
!!---------------------------------------------------------------------
!!     COMPUTE THE ADVECTIVE TERM: MERIDIONAL COMPONENT
!!---------------------------------------------------------------------
!!---------------------------------------------------------------------
!!     COMPUTE THE ADVECTIVE TERM: VERTICAL COMPONENT
!!---------------------------------------------------------------------
! 
!---------------------------------------------------------------------
!     COMPUTE THE VERTICAL VISCOSITY
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,K)= C0F * SQRT (UP (I,J,K)* UP (I,J,K) + VP (I, &
                            J,K)* VP (I,J,K))
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k,J,I,diff_u1,diff_v1,diff_u2,diff_v2,rkv,riv1,rkv1)
      DO K = 1,KM
      DO J = JSM,JEM
         DO I = 2,IMM
!lhl1204
#ifdef CANUTO
!       print*,akmU(i,j,k)
            if (k==1) then
               diff_v1 = SV (I,J)* OD0*(1-AIDIF)
               diff_u1 = SU (I,J)* OD0*(1-AIDIF)
            else
               diff_v1= AKMU(I,J,K-1)*(1-AIDIF)*(VP(I,J,K-1)- VP(I,J,K))*ODZT(K)*VIV(I,J,K)+ &
                        (1.0D0- VIV (I,J,K))* WKA (I,J,K -1)*(1-AIDIF) &
                      * (-SNLAT(J)*UP(I,J,K-1)*SAG+VP(I,J,K-1)*CAG)
               diff_u1= AKMU(I,J,K-1)*(1-AIDIF)* (UP(I,J,K-1)-UP(I,J,K))* ODZT (K)*VIV (I,J,K) + &
                       (1.0D0- VIV (I,J,K))* WKA (I,J,K -1)*(1-AIDIF) &
                      *(UP(I,J,K-1)*CAG+SNLAT(J)*VP(I,J,K-1)* SAG)
            end if
            if (k==km) then
               diff_v2= WKA (I,J,KM)* ( - SNLAT (J)* UP (I,J,KM)        &
                        * SAG + VP (I,J,KM)* CAG)*(1-AIDIF)
               diff_u2= WKA (I,J,KM)* ( UP (I,J,KM)* CAG + SNLAT (J)    &
                        * VP (I,J,KM)* SAG)*(1-AIDIF)
            else
               diff_v2= AKMU(I,J,K)*(1-AIDIF)*(VP(I,J,K)- VP(I,J,K+1))*ODZT(K+1)*VIV(I,J,K+1)+ &
                        (1.0D0- VIV (I,J,K+1))* WKA (I,J,K)*(1-AIDIF) &
                      * (-SNLAT(J)*UP(I,J,K)*SAG+VP(I,J,K)*CAG)
               diff_u2= AKMU(I,J,K)*(1-AIDIF)*(UP(I,J,K)-UP(I,J,K+1))* ODZT(K+1)*VIV (I,J,K+1) + &
                       (1.0D0- VIV (I,J,K+1))* WKA (I,J,K)*(1-AIDIF) &
                      *(UP(I,J,K)*CAG+SNLAT(J)*VP(I,J,K)* SAG)
            end if
#else
            RKV = AMV
            RKV1= AMV
#ifdef SPMD
            IF (J_global(j) >= RUST.AND.J_global(j) <= RUEND.AND.K>=2)THEN
#else
            IF (J >= RUST.AND.J <= RUEND.AND.K>=2)THEN
#endif
!        depended on Richardson number
            IF (riu (i,j,k -1) < 0.0)THEN
                RKV = visc_cbu_limit
            ELSE
                riv1 = 1.0D0/ (1.0D0+5.0D0* riu (i,j,k -1))
                RKV = fricmx * riv1* riv1+ visc_cbu_back
            END IF
                IF (k == 2.AND.RKV < wndmix) RKV = wndmix
            END IF
!
#ifdef SPMD
            IF (J_global(j) >= RUST.AND.J_global(j) <= RUEND)THEN
#else
            IF (J >= RUST.AND.J <= RUEND)THEN
#endif
            IF (riu (i,j,k) < 0.0)THEN
                RKV1= visc_cbu_limit
            ELSE
                riv1 = 1.0D0/ (1.0D0+5.0D0* riu (i,j,k))
                RKV1= fricmx * riv1* riv1+ visc_cbu_back
            END IF
                IF (k == 1.AND.RKV1< wndmix) RKV1= wndmix
            END IF
!
            AKMU(I,J,K)=RKV1
!
            if (k==1) then
               diff_v1 = SV (I,J)* OD0
               diff_u1 = SU (I,J)* OD0
            else
               diff_v1= RKV*(VP(I,J,K-1)- VP(I,J,K))*ODZT(K)*VIV(I,J,K)+ &
                        (1.0D0- VIV (I,J,K))* WKA (I,J,K -1) &
                      * (-SNLAT(J)*UP(I,J,K-1)*SAG+VP(I,J,K-1)*CAG)
               diff_u1= RKV * (UP(I,J,K-1)-UP(I,J,K))* ODZT (K)*VIV (I,J,K) + &
                       (1.0D0- VIV (I,J,K))* WKA (I,J,K -1) &
                      *(UP(I,J,K-1)*CAG+SNLAT(J)*VP(I,J,K-1)* SAG)   
            end if
            if (k==km) then
               diff_v2= WKA (I,J,KM)* ( - SNLAT (J)* UP (I,J,KM)        &
                        * SAG + VP (I,J,KM)* CAG)
               diff_u2= WKA (I,J,KM)* ( UP (I,J,KM)* CAG + SNLAT (J)    &
                        * VP (I,J,KM)* SAG)
            else
               diff_v2= RKV1*(VP(I,J,K)- VP(I,J,K+1))*ODZT(K+1)*VIV(I,J,K+1)+ &
                        (1.0- VIV (I,J,K+1))* WKA (I,J,K) &
                      * (-SNLAT(J)*UP(I,J,K)*SAG+VP(I,J,K)*CAG)
               diff_u2= RKV1* (UP(I,J,K)-UP(I,J,K+1))* ODZT(K+1)*VIV (I,J,K+1) + &
                       (1.0- VIV (I,J,K+1))* WKA (I,J,K) &
                      *(UP(I,J,K)*CAG+SNLAT(J)*VP(I,J,K)* SAG)   
            end if
!
!
#endif
!lhl1204
            DLV (I,J,K) = DLV (I,J,K) + ODZP (K)* (diff_v1-diff_v2)
            DLU (I,J,K) = DLU (I,J,K) + ODZP (K)* (diff_u1-diff_u2)
        END DO
      END DO
      END DO

#ifdef SPMD
      call exch_boundary(akmu(1,1,1),km)
#endif
      if (nx_proc.ne.1) then
      n2=mod(mytid,nx_proc)
      do k=1,km
      do j=1,jmt
        if (n2==0) then
         akmu(1,j,k)=0.0
        end if
        if (n2==nx_proc-1)then
         akmu(imt,j,k)=0.0
        end if
      end do
      end do
      end if

!
!lhl0711#ifdef SPMD
!lhl0711      call exch_boundary(akmu,km)
!lhl0711#endif
      deallocate(riu) 
 
 
!---------------------------------------------------------------------
!     COMPUTE THE HORIZONTAL VISCOSITY
!---------------------------------------------------------------------
#if ( defined SMAG)
 
      DO K = 1,KM
 
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = 1,JMT
            DO I = 1,IMT
               WKA (I,J,11)= UP (I,J,K)
               WKA (I,J,12)= VP (I,J,K)
            END DO
         END DO
!
         CALL SMAG2 (K)
!
#if (defined SMAG_FZ )
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,16)= (0.5D0* OUX(J)*(WKA(I +1,J,7)-WKA(I-1,J,7)) &
               -0.5D0* R2E (J)* WKA (I,J,8) +0.5D0* R2F (J)* WKA (I,J -1,8) &
                            -0.5D0* R2E (J +1)* WKA (I,J +1,8) +0.5D0* R2F (&
                             J +1)* WKA (I,J,8))* VIV (I,J,K)
               WKA (I,J,17)= (0.5D0* OUX (J)*(WKA(I+1,J,9)-WKA(I-1,J,9)) &
               -0.5D0* R3E (J)* WKA (I,J,10) +0.5D0* R3F (J)* WKA (I,J -1,10) &
               -0.5D0* R3E (J +1)* WKA (I,J +1,10) +0.5D0* R3F (J +1)* WKA (&
                           I,J,10) &
               + R4E (J)* WKA (I,J,7))* VIV (I,J,K)         
            END DO
         END DO
 
#else

!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,16)= (0.5D0* OUX (J)* (WKA (I +1,J,7) - WKA (I -1,&
                           J,7)) &
               - R2E (J)* WKA (I,J +1,8) + R2F (J)* WKA (I,J -1,8))     &
                * VIV (I,J,K)  
               WKA (I,J,17)= (0.5D0* OUX (J)* (WKA (I +1,J,9) - WKA (I -1,&
                           J,9)) &
               - R3E (J)* WKA (I,J +1,10) + R3F (J)* WKA (I,J -1,10) &
               + R4E (J)* WKA (I,J,7))* VIV (I,J,K)   
            END DO
         END DO
 
#endif
 
 
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               DLU (I,J,K)= DLU (I,J,K) + WKA (I,J,16)
               DLV (I,J,K)= DLV (I,J,K) + WKA (I,J,17)
            END DO
         END DO
      END DO
!
!!!!!!!
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               DLV (I,J,K)= DLV (I,J,K) + 0.2D0*AM (J)* &
               (R1D (J)* (VP (I,J +1,K) - VP (I,J,K)) &
               - R1C (J)* (VP (I,J,K) - VP (I,J -1,K)) &
               + SOUX (J)* (VP (I +1,J,K) -2.0D0* VP (I,J,K) + VP (I -1,J,K)) &
               + CV1 (J)* VP (I,J,K) &
               - CV2 (J)* (UP (I +1,J,K) - UP (I -1,J,K)))     
               DLU (I,J,K)= DLU (I,J,K) + 0.2D0*AM (J)* &
               (R1D (J)* (UP (I,J +1,K) - UP (I,J,K)) &
               - R1C (J)* (UP (I,J,K) - UP (I,J -1,K)) &
               + SOUX (J)* (UP (I +1,J,K) -2.0D0* UP (I,J,K) + UP (I -1,J,K)) &
               + CV1 (J)* UP (I,J,K) &
               + CV2 (J)* (VP (I +1,J,K) - VP (I -1,J,K)))    
            END DO
         END DO
       END DO
!!!!!!!
#else
      
#if (defined BIHAR)

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               TMP1 (I,J,K)= AM3 (I,J,K)* &
               (R1D (J)* (VP (I,J +1,K) - VP (I,J,K)) &
               - R1C (J)* (VP (I,J,K) - VP (I,J -1,K)) &
               + SOUX (J)* (VP (I +1,J,K) -2.0D0* VP (I,J,K) + VP (I -1,J,K)) &
               + CV1 (J)* VP (I,J,K) &
               - CV2 (J)* (UP (I +1,J,K) - UP (I -1,J,K)))* VIV (I,J,K) 
               TMP2 (I,J,K)= AM3 (I,J,K)* &
               (R1D (J)* (UP (I,J +1,K) - UP (I,J,K)) &
               - R1C (J)* (UP (I,J,K) - UP (I,J -1,K)) &
               + SOUX (J)* (UP (I +1,J,K) -2.0D0* UP (I,J,K) + UP (I -1,J,K)) &
               + CV1 (J)* UP (I,J,K) &
               + CV2 (J)* (VP (I +1,J,K) - VP (I -1,J,K)))* VIV (I,J,K)
            END DO
         END DO
      END DO
      if (nx_proc == 1) then
!$OMP PARALLEL DO PRIVATE (K,J)
         do k=1,km
         do j=jsm,jem
            TMP1 (1,J,K)= TMP1 (IMM,J,K)
            TMP1 (IMT,J,K)= TMP1 (2,J,K)
            TMP2 (1,J,K)= TMP2 (IMM,J,K)
            TMP2 (IMT,J,K)= TMP2 (2,J,K)
        end do
        end do
      end if
!YU
#ifdef SPMD
      call exch_boundary(tmp1,km)
      call exch_boundary(tmp2,km)
#endif 
!YU
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               DLV (I,J,K)= DLV (I,J,K) + &
               (R1D (J)* (TMP1 (I,J +1,K) - TMP1 (I,J,K)) &
               - R1C (J)* (TMP1 (I,J,K) - TMP1 (I,J -1,K)) &
               + SOUX (J)* (TMP1 (I +1,J,K) -2.0D0* TMP1 (I,J,K) + TMP1 ( &
                           I -1,J,K))&
               + CV1 (J)* TMP1 (I,J,K) &
               - CV2 (J)* (TMP2 (I +1,J,K) - TMP2 (I -1,J,K)))* VIV (I,J,K)
               DLU (I,J,K)= DLU (I,J,K) + &
               (R1D (J)* (TMP2 (I,J +1,K) - TMP2 (I,J,K)) &
               - R1C (J)* (TMP2 (I,J,K) - TMP2 (I,J -1,K)) &
               + SOUX (J)* (TMP2 (I +1,J,K) -2.0D0* TMP2 (I,J,K) + TMP2 ( &
                           I -1,J,K))&
               + CV1 (J)* TMP2 (I,J,K) &
               + CV2 (J)* (TMP1 (I +1,J,K) - TMP1 (I -1,J,K)))* VIV (I,J,K) 
            END DO
         END DO
      END DO
      deallocate(tmp1,tmp2)
#else
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               DLV (I,J,K)= DLV (I,J,K) + AM3 (I,J,K)* &
               (R1D (J)* (VP (I,J +1,K) - VP (I,J,K)) &
               - R1C (J)* (VP (I,J,K) - VP (I,J -1,K)) &
               + SOUX (J)* (VP (I +1,J,K) -2.0D0* VP (I,J,K) + VP (I -1,J,K)) &
               + CV1 (J)* VP (I,J,K) &
               - CV2 (J)* (UP (I +1,J,K) - UP (I -1,J,K)))     
               DLU (I,J,K)= DLU (I,J,K) + AM3 (I,J,K)* &
               (R1D (J)* (UP (I,J +1,K) - UP (I,J,K)) &
               - R1C (J)* (UP (I,J,K) - UP (I,J -1,K)) &
               + SOUX (J)* (UP (I +1,J,K) -2.0D0* UP (I,J,K) + UP (I -1,J,K)) &
               + CV1 (J)* UP (I,J,K) &
               + CV2 (J)* (VP (I +1,J,K) - VP (I -1,J,K)))    
            END DO
         END DO
      END DO
 
#endif
#endif
 
!---------------------------------------------------------------------
!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!---------------------------------------------------------------------
 
      allocate(dlub(imt,jmt),dlvb(imt,jmt),stat=ierr)
      CALL VINTEG (DLU,DLUB)
      CALL VINTEG (DLV,DLVB)
!
     if ( nx_proc == 1 ) then 
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
 

 
 
!---------------------------------------------------------------------
!     VERTICAL INTEGRATION
!---------------------------------------------------------------------
 
!Yu   if (ierr /= 0) then
!Yu      write(6,*)'allocation error---tmp1,tmp2'
!Yu      stop
!Yu   end if


      RETURN
      END SUBROUTINE READYC
 
 
