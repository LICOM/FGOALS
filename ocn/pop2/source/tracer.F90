!  CVS: $Id: tracer.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE TRACER
!     =================
#include <def-undef.h>
use precision_mod 
use param_mod
use pconst_mod
use tracer_mod
use work_mod
use dyn_mod
use isopyc_mod
use forc_mod
use pmix_mod
#ifdef SPMD
use msg_mod
#endif
 
      IMPLICIT NONE
 
      integer     :: n2
      REAL(r8)    :: AIDIF,C2DTTS,AA,FAW,FIW,ALF,RNCC,ABC,fil_lat1,fil_lat2

!Xiao Chan (Hereinafter XC for short)
#if (defined TSPAS)
      real(r8)    :: LAMDA,wt1,wt2,adv_z
      real(r8),dimension(:,:,:), allocatable :: adv_xy1,adv_xy2,adv_xy3,adv_xy4
!     real(r8),dimension(:,:,:), allocatable :: adv_xy1,adv_xy2,adv_xy3,adv_xy4, &
!               adv_x0,adv_y0,adv_c1,adv_c2,atmax,atmin,at0,uaa,vaa,adv_xx,adv_yy
      real(r8),dimension(imt,jmt,km) :: uaa,vaa, &
                adv_x0,adv_y0,adv_c1,adv_c2,atmax,atmin,at0,adv_xx,adv_yy
      real(r8),dimension(:,:,:) , allocatable :: adv_zz,atz, adv_za,adv_zb1,adv_zb2,adv_zc,atmaxz,atminz
#else
      real(r8)    :: LAMDA,wt1,wt2,adv_y,adv_x,adv_z,adv_x1,adv_x2
#endif
!XC
 
!---------------------------------------------------------------------
!     SET LOCAL CONSTANT
!---------------------------------------------------------------------
      deallocate(dlu,dlv,gg,ric,rict)
 
      allocate(stf(imt,jmt),tf(imt,jmt,km))
      allocate(wkb(imt,jmt,km),wkc(imt,jmt,km),wkd(imt,jmt,km))
!lhl      AIDIF = 0.5*FLOAT(ISOP)
!
!---------------------------------------------------------------------
!      Define the threthold latitute for zonal smoother
       fil_lat1=63.0D0
       fil_lat2=63.0D0

#if (defined ISO)
      AIDIF = 0.5D0
#else
      AIDIF = 0.0D0
#endif
 
      LAMDA=1.0D0/(15.D0*86400.D0)
      RNCC = 1.0D0/ FLOAT (NCC)
      IF (IST >= 1)THEN
         C2DTTS = DTS *2.0D0
         AA = 0.5D0
      ELSE
         C2DTTS = DTS
         AA = 0.0D0
      END IF
 

!---------------------------------------------------------------------
!     PREPARATION FOR VERTICAL ADVECTIVE TERM
!---------------------------------------------------------------------
#ifdef COUP
!LPF 20120823
#ifdef BACKMX
!M
!$OMP PARALLEL DO PRIVATE (k,J,I)
      do k=1,22
      do j=1,jmt
      do i=1,imt
         if (basin(i,j_global(j)) == 1 .or. basin(i,j_global(j)) == 2) then
            if ( ahv_back(j_global(j)) > 1.0D-5 .and. at(i,j,1,1) < 4.0 ) then
               akt(i,j,k,1) = max ( akt(i,j,k,1), ahv_back(j_global(j)))
               akt(i,j,k,2) = max ( akt(i,j,k,2), ahv_back(j_global(j)))
            end if
         end if
!
         if (basin(i,j_global(j)) == 4) then
            if ( ahv_back(j_global(j)) > 1.0D-5 .and. at(i,j,1,1) < 2.0 .and. k < 16) then
               akt(i,j,k,1) = max ( akt(i,j,k,1), ahv_back(j_global(j)))
               akt(i,j,k,2) = max ( akt(i,j,k,2), ahv_back(j_global(j)))
            end if
         end if
      end do
      end do
      end do
#endif
#endif 
!LPF 20120823
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            H0F (I,J)= H0F (I,J)* ONBC
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            STF (I,J)= AA * H0F (I,J) + (1.0D0- AA)* H0L (I,J)
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UTF (I,J,K)= UTF (I,J,K)* ONCC
               VTF (I,J,K)= VTF (I,J,K)* ONCC
            END DO
         END DO
      END DO

#ifdef DEBUG
      call chk_var2d(utf,"ut",0)
      call chk_var2d(vtf,"vt",0)
#endif

 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKD (I,J,K)= AA * UTF (I,J,K) + (1.0D0- AA)* UTL (I,J,K)
               WKC (I,J,K)= AA * VTF (I,J,K) + (1.0D0- AA)* VTL (I,J,K)
#if ( defined SMAG)
               UTL (I,J,K)= AA * UTF (I,J,K) + (1.0D0- AA)* UTL (I,J,K)
               VTL (I,J,K)= AA * VTF (I,J,K) + (1.0D0- AA)* VTL (I,J,K)
#endif
            END DO
         END DO
      END DO

#ifdef DEBUG
      call chk_var2d(wkd,"kd",0)
      call chk_var2d(wkc,"kc",0)
#endif
!XC
#if (defined TSPAS)
!     allocate (uaa(imt,jmt,km),vaa(imt,jmt,km)) 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               uaa(i,j,k)=0.5D0*WKD (i,j,k) + 0.5D0*WKD (i,j-1,k)
               vaa(i,j,k)=0.5D0*WKC (i,j,k) + 0.5D0*WKC (i+1,j,k)
            END DO
         END DO
      END DO
         if (nx_proc==1) then
!$OMP PARALLEL DO PRIVATE (K,J)
         DO K=1,KM
         DO J = JST,JET
               uaa (1,J,K) = uaa(imm,J,K)
               uaa (imt,J,K) = uaa(2,J,K)
               vaa (1,J,K) = vaa(imm,J,K)
               vaa (imt,J,K) = vaa(2,J,K)
         END DO
         END DO
         end if

#ifdef SPMD
       call exch_boundary(uaa(1,1,1),km)
       call exch_boundary(vaa(1,1,1),km)
#endif
#endif
!XC
 
      CALL UPWELL (WKD,WKC,STF)
#ifdef DEBUG
      call chk_var2d(ws,"ws",1)
#endif
 
#if (defined NODIAG)
 
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
 
!---------------------------------------------------------------------
!     PREPARATION FOR HORIZONAL ADVECTIVE TERM
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 1,IMT
               UTL (I,J,K)= 0.25D0* OTX (J)* (WKD (I,J,K) + WKD (I,J -1,K))
            END DO
         END DO
      END DO
#ifdef DEBUG
      call chk_var2d(utl,"tl",1)
#endif
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               WKD (I,J,K)= R2A (J)* (WKC (I,J,K) + WKC (I +1,J,K))
               WKB (I,J,K)= R2B (J)* (WKC (I,J -1,K) + WKC (I +1,J -1,K))
            END DO
         END DO
      END DO
#ifdef DEBUG
      call chk_var2d(wkd,"kd",1)
      call chk_var2d(wkd,"kc",1)
#endif
 
!-----------------------------------------------------------------------
!     PREPARATION FOR ISOPYCNAL DIFFUSION & ADVECTION
!-----------------------------------------------------------------------
#if (defined ISO)
      CALL ISOPYC
#endif
 
!@@@  COMPUTING DIFFUSION COEFFICIENT
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
#if (defined ISO)
               WKC (I,J,K) = AHV + AHISOP * K3 (I,K,J,3)
#else
               WKC (I,J,K) = AHV
#endif
!
#if (!defined CANUTO)
            AKT(I,J,K,1)=WKC(I,J,K)
            AKT(I,J,K,2)=WKC(I,J,K)
#endif
!
            END DO
         END DO
      END DO
 
!     AT LOW LATITUDE, DIFFUSION DEPEND ON RICHARDSON NUMBER
 
!$OMP PARALLEL DO PRIVATE (K,J,I,ABC,ALF)
      DO K = 1,KMM1
      do j=  jsm,jem
#ifdef SPMD
         if (j_global(j)>=rtst.and.j_global(j)<=rtend) then
#else
         if (j>=rtst.and.j<=rtend) then
#endif
            DO I = 1,IMT
               IF (rit (i,j,k) < 0.0) THEN
                  ABC = diff_cbt_limit
               ELSE
                  ALF = 1.0D0/ (1.0D0+5.0D0* rit (i,j,k)* RNCC)
                  ABC = fricmx * ALF **3+ diff_cbt_back
               END IF
               IF (k == 1.AND.ABC < wndmix) ABC = wndmix
!lhl        WKC(I,J,K) = ABC+AHISOP*K3(I,K,J,3)*FLOAT(ISOP)
#if (defined ISO)
               WKC (I,J,K) = ABC + AHISOP * K3 (I,K,J,3)
#else
               WKC (I,J,K) = ABC
#endif
!
#if (!defined CANUTO)
            AKT(I,J,K,1)=WKC(I,J,K)
            AKT(I,J,K,2)=WKC(I,J,K)
#endif
!
            END DO
         end if
      end do
      end do
!
#if (!defined CANUTO)
#ifdef SPMD
      call exch_boundary(akt(1,1,1,1),km)
      call exch_boundary(akt(1,1,1,2),km)
#endif
#endif
!
#ifdef DEBUG
      call chk_var2d(wkc,"kc",1)
#endif
!-----------------------------------------------------------------------
!     SOLVE FOR ONE TRACER AT A TIME
!-----------------------------------------------------------------------
!     NTRA = 1 => TEMPERATURE
!     NTRA = 2 => SALINITY
 
      DO N = 1,NTRA

#ifdef CANUTO      
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = JST,JET          
            DO I = 2,IMM
#if (defined ISO)       
               WKC (I,J,K) = AKT(I,J,K,N) + AHISOP * K3 (I,K,J,3) 
#else            
               WKC (I,J,K) = AKT(I,J,K,N) 
#endif         
            END DO
         END DO
      END DO
!
         if (nx_proc==1) then
!$OMP PARALLEL DO PRIVATE (K,J)
         DO K=1,KM
         DO J = JST,JET          
               WKC (1,J,K) = WKC(imm,J,K)
               WKC (imt,J,K) = WKC(2,J,K)
         END DO
         END DO
         end if
#ifdef SPMD
      call exch_boundary(wkc,km)
#endif
   
              
#endif
!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM 
!---------------------------------------------------------------------

#if (defined TSPAS)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!XC!!$OMP PARALLEL DO PRIVATE (K,J,I,adv_x,adv_x1,adv_x2,adv_y,adv_z,wt1,wt2)
!XiaoChan(XC)
!Use Two-step Shape-Preserving Advection Scheme (TSPAS)

!     allocate ( adv_x0(imt,jmt,km),adv_y0(imt,jmt,km),adv_c1(imt,jmt,km), adv_c2(imt,jmt,km), adv_xx(imt,jmt,km), adv_yy(imt,jmt,km))
      allocate ( adv_xy1(imt,jmt,km), adv_xy2(imt,jmt,km), adv_xy3(imt,jmt,km), adv_xy4(imt,jmt,km))
!     allocate ( at0(imt,jmt,km) )
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,km
        DO J = JSM,JEM
          DO I = 2,IMM
            adv_x0(i,j,k)=0.5*otx (j)*(at(i+1,j,k,n)+at(i,j,k,n))*uaa(i+1,j,k) &
                       -0.5*otx (j)*(at(i,j,k,n)+at(i-1,j,k,n))*uaa(i,j,k)
            adv_y0(i,j,k)=0.5*(at(i,j+1,k,n)+at(i,j,k,n))*vaa(i,j,k)*4*R2A(j) &
                       -0.5*(at(i,j,k,n)+at(i,j-1,k,n))*vaa(i,j-1,k)*4*R2B(j)
            adv_xy1(i,j,k)=-0.5*otx(j)*dtdx(j)*(at(i+1,j,k,n)-at(i,j,k,n))*  &
                                  uaa(i+1,j,k)*uaa(i+1,j,k)
            adv_xy2(i,j,k)= 0.5*otx(j)*dtdx(j)*(at(i,j,k,n)-at(i-1,j,k,n))*  &
                                  uaa(i,j,k)*uaa(i,j,k)
            adv_xy3(i,j,k)=-0.5*dtdy(j)*(at(i,j+1,k,n)-at(i,j,k,n))*  &
                                  vaa(i,j,k)*vaa(i,j,k)*4*R2A(j)*RAA(j)
            adv_xy4(i,j,k)= 0.5*dtdy(j-1)*(at(i,j,k,n)-at(i,j-1,k,n))*  &
                                  vaa(i,j-1,k)*vaa(i,j-1,k)*4*R2B(j)*RBB(j)
            adv_c1(i,j,k)=-AT(i,j,k,n)*otx(j)*(UAA(i+1,j,k)-UAA(i,j,k))
            adv_c2(i,j,k)=-AT(i,j,k,n)*4*(R2A(j)*VAA(i,j,k)-R2B(j)*VAA(i,j-1,k))

            adv_xx(i,j,k)=-(adv_x0(i,j,k)+adv_xy1(i,j,k)+adv_xy2(i,j,k)+adv_c1(i,j,k))
            adv_yy(i,j,k)=-(adv_y0(i,j,k)+adv_xy3(i,j,k)+adv_xy4(i,j,k)+adv_c2(i,j,k))

            at0(i,j,k)=at(i,j,k,n)+(adv_xx(i,j,k)+adv_yy(i,j,k))*dts
          ENDDO
        ENDDO
      ENDDO

!Circle condition
      if (nx_proc==1) then
!$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1,km
        DO J = JSM,JEM
           at0(1  ,j,k)=at0(IMM,j,k)
           at0(IMT,j,k)=at0(2  ,j,k)
        ENDDO
      ENDDO
      end if
#ifdef SPMD
       call exch_boundary(at0,km)
#endif

!XC !!$OMP PARALLEL DO PRIVATE (I)
!        DO I = 1,JMT
!           at0(I,jsm-1)=at0(I,jsm)
!           at0(I,jem+1)=at0(I,jem)
!        ENDDO
!Compute Max and Min
!     allocate ( atmax(imt,jmt,km), atmin(imt,jmt,km) )
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,km
        DO J = JSM,JEM
          DO I = 2,IMM
            atmax(i,j,k)=max(at(i,j,k,n),at(i,j-1,k,n),at(i,j+1,k,n), &
                           at(i-1,j,k,n),at(i+1,j,k,n) )
            atmin(i,j,k)=min(at(i,j,k,n),at(i,j-1,k,n),at(i,j+1,k,n), &
                           at(i-1,j,k,n),at(i+1,j,k,n) )
          ENDDO
        ENDDO
      ENDDO

     if (nx_proc==1) then
!$OMP PARALLEL DO PRIVATE (K,J)
      DO K = 1,km
        DO J = JSM,JEM
           atmax(1  ,j,k)=atmax(IMM,j,k)
           atmin(IMT,j,k)=atmin(2  ,j,k)
        ENDDO
      ENDDO
     end if
!XC !!$OMP PARALLEL DO PRIVATE (K,I)
!      DO K = 1,km
!        DO I = 1,JMT
!           atmax(I,jsm-1,k)=atmax(I,jsm,k)
!           atmin(I,jem+1,k)=atmin(I,jem,k)
!        ENDDO
!      ENDDO

#ifdef SPMD
       call exch_boundary(atmax,km)
       call exch_boundary(atmin,km)
#endif

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,km
        DO J = JSM,JEM
          DO I = 2,IMM
            if (at0(i,j,k)>atmax(i,j,k).or.at0(i,j,k)<atmin(i,j,k)) then
              adv_xy1(i,j,k)=-0.5*otx(j)*(at(i+1,j,k,n)-at(i,j,k,n))* &
                                 abs(uaa(i+1,j,k))
              adv_xy2(i,j,k)= 0.5*otx(j)*(at(i,j,k,n)-at(i-1,j,k,n))* &
                                 abs(uaa(i,j,k))
              adv_xy3(i,j,k)=-0.5*(at(i,j+1,k,n)-at(i,j,k,n))* &
                                 abs(vaa(i,j,k))*4*R2A(j)
              adv_xy4(i,j,k)= 0.5*(at(i,j,k,n)-at(i,j-1,k,n))* &
                                 abs(vaa(i,j-1,k))*4*R2B(j)
            else
              if (at0(i+1,j,k)>atmax(i+1,j,k).or.at0(i+1,j,k)<atmin(i+1,j,k)) then
                adv_xy1(i,j,k)=-0.5*otx(j)*(at(i+1,j,k,n)-at(i,j,k,n))* &
                                 abs(uaa(i+1,j,k))
              endif
              if (at0(i-1,j,k)>atmax(i-1,j,k).or.at0(i-1,j,k)<atmin(i-1,j,k)) then
                adv_xy2(i,j,k)= 0.5*otx(j)*(at(i,j,k,n)-at(i-1,j,k,n))* &
                                 abs(uaa(i,j,k))
              endif
              if (at0(i,j+1,k)>atmax(i,j+1,k).or.at0(i,j+1,k)<atmin(i,j+1,k)) then
                adv_xy3(i,j,k)=-0.5*(at(i,j+1,k,n)-at(i,j,k,n))* &
                                 abs(vaa(i,j,k))*4*R2A(j)
              endif
              if (at0(i,j-1,k)>atmax(i,j-1,k).or.at0(i,j-1,k)<atmin(i,j-1,k)) then
                adv_xy4(i,j,k)= 0.5*(at(i,j,k,n)-at(i,j-1,k,n))* &
                                 abs(vaa(i,j-1,k))*4*R2B(j)
              endif
            endif

            adv_xx(i,j,k)=-(adv_x0(i,j,k)+adv_xy1(i,j,k)+adv_xy2(i,j,k)+adv_c1(i,j,k))
            adv_yy(i,j,k)=-(adv_y0(i,j,k)+adv_xy3(i,j,k)+adv_xy4(i,j,k)+adv_c2(i,j,k))

                  TF (I,J,K)= adv_xx(i,j,k)+adv_yy(i,j,k)

                  ax(i,j,k,N) = adv_xx(i,j,k)
                  ay(i,j,k,N) = adv_yy(i,j,k)

          END DO
        END DO
      END DO

      deallocate ( adv_xy1,adv_xy2,adv_xy3,adv_xy4)
!     deallocate ( atmax, atmin, uaa, vaa, adv_xy1,adv_xy2,adv_xy3,adv_xy4)
!     deallocate ( adv_x0,adv_y0,adv_c1, adv_c2, adv_xx, adv_yy, at0)
#ifdef SPMD
       call exch_boundary(ax(1,1,1,n),km)
       call exch_boundary(ay(1,1,1,n),km)
#endif

!Use Shape-Preserving Scheme in vertial direction
      allocate ( adv_zz(imt,jmt,km), adv_za(imt,jmt,km),adv_zb1(imt,jmt,km),adv_zb2(imt,jmt,km), adv_zc(imt,jmt,km), atmaxz(imt,jmt,km), atminz(imt,jmt,km), atz(imt,jmt,km))
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,km
        DO J = JSM,JEM
          DO I = 2,IMM
            if (k==1) then
              adv_za (i,j,k)=-0.5*ODZP(1)*WS(I,J,2)*(AT(I,J,2,N)+AT(I,J,1,N))
              adv_zb1(i,j,k)=0
              adv_zb2(i,j,k)= 0.5*ODZP(1)*WS(I,J,2)*WS(I,J,2)*ODZT(2) &
                                *(AT(I,J,1,N)-AT(I,J,2,N))
              adv_zc (i,j,k)=     ODZP(1)*AT(I,J,1,N)*WS(I,J,2)
              atmaxz(i,j,k)=max(at(i,j,1,n),at(i,j,2,n))
              atminz(i,j,k)=min(at(i,j,1,n),at(i,j,2,n))
              atz(i,j,k)=at(i,j,k,n)-(adv_za(i,j,k)+adv_zb1(i,j,k) &
                                     +adv_zb2(i,j,k)+adv_zc(i,j,k))*dts
            elseif (k==km) then
              adv_za (i,j,k)= 0.5*ODZP(km)*WS(I,J,km)*(AT(I,J,km,N)+AT(I,J,km-1,N))
              adv_zb1(i,j,k)=-0.5*ODZP(km)*WS(I,J,km)*WS(I,J,km)*ODZT(km  ) &
                                *(AT(I,J,km-1,N)-AT(I,J,km,N))
              adv_zb2(i,j,k)=0
              adv_zc (i,j,k)=    -ODZP(km)*AT(I,J,km,N)*WS(I,J,km)
              atmaxz(i,j,k)=max(at(i,j,km-1,n),at(i,j,km,n))
              atminz(i,j,k)=min(at(i,j,km-1,n),at(i,j,km,n))
              atz(i,j,k)=at(i,j,k,n)-(adv_za(i,j,k)+adv_zb1(i,j,k) &
                                     +adv_zb2(i,j,k)+adv_zc(i,j,k))*dts
            else
              adv_za (i,j,k)= 0.5*ODZP(k)*WS(I,J,k  )*(AT(I,J,k,N)+AT(I,J,k-1,N)) &
                         -0.5*ODZP(k)*WS(I,J,k+1)*(AT(I,J,k,N)+AT(I,J,k+1,N))
              adv_zb1(i,j,k)=-0.5*ODZP(k)*WS(I,J,k  )*WS(I,J,k  )*ODZT(k  ) &
                                *(AT(I,J,k-1,N)-AT(I,J,k,N))
              adv_zb2(i,j,k)= 0.5*ODZP(k)*WS(I,J,k+1)*WS(I,J,k+1)*ODZT(k+1) &
                                *(AT(I,J,k,N)-AT(I,J,k+1,N))
              adv_zc (i,j,k)=    -ODZP(k)*AT(I,J,k,N)*(WS(I,J,k)-WS(I,J,k+1))
              atmaxz (i,j,k)=max(at(i,j,k-1,n),at(i,j,k,n),at(i,j,k+1,n))
              atminz (i,j,k)=min(at(i,j,k-1,n),at(i,j,k,n),at(i,j,k+1,n))
              atz(i,j,k)=at(i,j,k,n)-(adv_za(i,j,k)+adv_zb1(i,j,k) &
                                     +adv_zb2(i,j,k)+adv_zc(i,j,k))*dts
            endif
          END DO
        END DO
      END DO

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,km
        DO J = JSM,JEM
          DO I = 2,IMM
            if (k==1) then
              if(atz(i,j,k+1)>atmaxz(i,j,k+1).or.atz(i,j,k+1)<atminz(i,j,k+1).or. &
                 atz(i,j,k)>atmaxz(i,j,k).or.atz(i,j,k)<atminz(i,j,k)) then
                adv_zb2(i,j,k)= 0.5*abs(WS(I,J,k+1))*ODZT(k+1) &
                                 *(AT(I,J,k,N)-AT(I,J,k+1,N))
              endif
             elseif (k==km) then
               if(atz(i,j,k-1)>atmaxz(i,j,k-1).or.atz(i,j,k-1)<atminz(i,j,k-1).or. &
                  atz(i,j,k  )>atmaxz(i,j,k  ).or.atz(i,j,k  )<atminz(i,j,k  )) then
                 adv_zb1(i,j,k)=-0.5*abs(WS(I,J,k  ))*ODZT(k  ) &
                                  *(AT(I,J,k-1,N)-AT(I,J,k,N))
               endif
             else
               if(atz(i,j,k)>atmaxz(i,j,k).or.atz(i,j,k)<atminz(i,j,k)) then
                 adv_zb1(i,j,k)=-0.5*abs(WS(I,J,k  ))*ODZT(k  ) &
                                  *(AT(I,J,k-1,N)-AT(I,J,k,N))
                 adv_zb2(i,j,k)= 0.5*abs(WS(I,J,k+1))*ODZT(k+1) &
                                  *(AT(I,J,k,N)-AT(I,J,k+1,N))
               else
                 if(atz(i,j,k+1)>atmaxz(i,j,k+1).or.atz(i,j,k+1)<atminz(i,j,k+1)) then
                   adv_zb2(i,j,k)= 0.5*abs(WS(I,J,k+1))*ODZT(k+1) &
                                  *(AT(I,J,k,N)-AT(I,J,k+1,N))
                 endif
                 if(atz(i,j,k-1)>atmaxz(i,j,k-1).or.atz(i,j,k-1)<atminz(i,j,k-1)) then
                   adv_zb1(i,j,k)=-0.5*abs(WS(I,J,k  ))*ODZT(k  ) &
                                  *(AT(I,J,k-1,N)-AT(I,J,k,N))
                 endif
               endif
             endif

             adv_zz(i,j,k)=-(adv_za(i,j,k)+adv_zb1(i,j,k)+adv_zb2(i,j,k)+adv_zc(i,j,k))
             atz(i,j,k)=at(i,j,k,n)+adv_zz(i,j,k)*dts

                  TF (I,J,K)= TF (I,J,K) + adv_zz(i,j,k)
!
                  az(i,j,k,N) = adv_zz(i,j,k)
!
          END DO
        END DO
      END DO

      deallocate ( adv_za,adv_zb1,adv_zb2, adv_zc, atmaxz, atminz, atz, adv_zz)
!      call chk_var(tf,km,trave)
#ifdef SPMD
       call exch_boundary(az(1,1,1,n),km)
#endif

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#else
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Use the centered-in-time centered-in-space advection scheme

!$OMP PARALLEL DO PRIVATE (K,J,I,adv_x,adv_x1,adv_x2,adv_y,adv_z,wt1,wt2)
         DO K = 2,km-1
            DO J = JSM,JEM
               DO I = 2,IMM
                  adv_x1=(AT (I  ,J,K,N) - AT (I-1,J,K,N))* UTL (I,J,K)
                  adv_x2=(AT (I+1,J,K,N) - AT (I  ,J,K,N))* UTL (I+1,J,K)
                  adv_x = - (adv_x2+adv_x1)
                  adv_y=-(WKD(I,J,K)*(AT(I,J+1,K,N)-AT(I,J,K,N))+WKB(I,J,K)*(AT(I,J,K,N)-AT(I,J-1,K,N)))
                  wt1= WS(I,J,K)*(AT(I,J,K-1,N)-AT(I,J,K,N))
                  wt2= WS(I,J,K+1)*(AT(I,J,K,N)-AT(I,J,K+1,N))
                  adv_z=-0.5D0*(wt1+wt2)*ODZP(K)
                  TF (I,J,K)= adv_x+adv_y+adv_z
!
                  ax(i,j,k,N) = adv_x
                  ay(i,j,k,N) = adv_y
                  az(i,j,k,N) = adv_z
!
               END DO
            END DO
         END DO

!$OMP PARALLEL DO PRIVATE (J,I,adv_x,adv_x1,adv_x2,adv_y,adv_z)
        DO J=JSM,JEM
           DO I = 2,IMM
              adv_x1=(AT (I,J,1,N) - AT (I-1,J,1,N))* UTL (I,J,1)
              adv_x2=(AT (I+1,J,1,N) - AT (I,J,1,N))* UTL (I+1,J,1)
              adv_x = - (adv_x2+adv_x1)
              adv_y=-(WKD(I,J,1)*(AT(I,J+1,1,N)-AT(I,J,1,N))+WKB(I,J,1)*(AT(I,J,1,N)-AT(I,J-1,1,N)))
              adv_z= -0.5D0*WS(I,J,2)*(AT(I,J,1,N)-AT(I,J,2,N))*odzp(1)
              TF (I,J,1)= adv_x+adv_y+adv_z
!
                  ax(i,j,1,N) = adv_x
                  ay(i,j,1,N) = adv_y
                  az(i,j,1,N) = adv_z
!
           END DO
        END DO
!$OMP PARALLEL DO PRIVATE (J,I,adv_x,adv_x1,adv_x2,adv_y,adv_z)
        DO J=JSM,JEM
           DO I = 2,IMM
              adv_x1=(AT (I  ,J,km,N) - AT (I-1,J,km,N))* UTL (I,J,km)
              adv_x2=(AT (I+1,J,km,N) - AT (I,J,km,N))* UTL (I+1,J,km)
              adv_x = - (adv_x2+adv_x1)
              adv_y=-(WKD(I,J,km)*(AT(I,J+1,km,N)-AT(I,J,km,N))+WKB(I,J,km)*(AT(I,J,km,N)-AT(I,J-1,km,N)))
              adv_z= -0.5D0*WS(I,J,km)*(AT(I,J,km-1,N)-AT(I,J,km,N))*odzp(km)
              TF (I,J,km)= adv_x+adv_y+adv_z
!
                  ax(i,j,km,N) = adv_x
                  ay(i,j,km,N) = adv_y
                  az(i,j,km,N) = adv_z
!
           END DO
        END DO

#ifdef SPMD
       call exch_boundary(ax(1,1,1,n),km)
       call exch_boundary(ay(1,1,1,n),km)
       call exch_boundary(az(1,1,1,n),km)
#endif

#ifdef DEBUG
       call chk_var2d(tf,"ad",1)
#endif
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#endif
!-----------------------------------------------------------------------
!     COMPUTE THE ISOPYCNAL/DIPYCNAL MIXING
!-----------------------------------------------------------------------
!     XZ AND YZ ISOPYCNAL DIFFUSIVE FLUX ARE SOLVED EXPLICITLY;
!     WHILE ZZ COMPONENT WILL BE SOLVED IMPLICITLY.
 
 
#if (defined ISO)
         CALL ISOFLUX (N)
#else 
 
#if ( defined SMAG)
         CALL SMAG3
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = JSM,JEM
               WKI (1)= 0.0D0
               DO I = 2,IMT
                  WKI (I) = 0.5D0* (AH3 (I,J,K) + AH3 (I -1,J,K))* (ATB ( &
                           I,J,K,N) - ATB (I -1,J,K,N))* VIT (I,J,K)* VIT (I -1,J,K)
               END DO
               DO I = 2,IMM
                  TF (I,J,K) = TF (I,J,K) + SOTX (J)*(WKI(I+1)-WKI(I))
!
!                dx(i,j,k,n)= SOTX (J)*(WKI(I+1)-WKI(I))
!
               END DO
            END DO
         END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I,wt1,wt2)
         DO K = 1,KM
         DO J = JSM,JEM
         DO I = 2,IMM
            wt1= 0.5D0*(AH3(I,J,K)+AH3(I,J-1,K))*(ATB(I,J,K,N)-ATB(I,J-1,K,N))*VIT(I,J,K)*VIT(I,J-1,K)
            wt2= 0.5D0*(AH3(I,J+1,K)+AH3(I,J,K))*(ATB(I,J+1,K,N)-ATB(I,J,K,N))*VIT(I,J+1,K)*VIT(I,J,K)
            TF (I,J,K) = TF (I,J,K) + (R2D (J)* wt2 - R2C(J)* wt1)
!
!           dy(i,j,k,n)=(R2D(J)*wt2-R2C(J)*wt1)
!          ddy(i,j,k,n)=(R2A(J)*wt2-R2B(J)*wt1)
!
         END DO
         END DO
         END DO
#ifdef SPMD
!      call exch_boundary(dx(1,1,1,n),km,1,1)
!      call exch_boundary(dy(1,1,1,n),km,1,1)
#endif
 
#else
 
!-----------------------------------------------------------------------
!     COMPUTE THE EDDY-DIFFUSION TERM :  ZONAL COMPONENT
!-----------------------------------------------------------------------
 
#if (defined BIHAR)
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,IMT
                  WKA (I,J,K) = 0.0D0
               END DO
            END DO
         END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = JSM,JEM
               WKI (1)= 0.0D0
               DO I = 2,IMT
                  WKI (I) = (ATB (I,J,K,N) - ATB (I -1,J,K,N))* VIT (I,J,K)* VIT (I -1,J,K)
               END DO
 
               DO I = 2,IMM
                  WKA (I,J,K) = AH3 (I,J,K)* SOTX (J)* (WKI (I +1) - WKI (I))
               END DO
            END DO
         END DO
 
        if (nx_proc==1) then
         DO K = 1,KM
             DO J = JSM,JEM
               WKA (1,J,K)= WKA (IMM,J,K)
               WKA (IMT,J,K)= WKA (2,J,K)
             END DO
        END DO
        end if
#ifdef SPMD
         call exch_boundary(WKA(1,1,1),km)
#endif

!!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = JSM,JEM
               WKI (1)= 0.0
               DO I = 2,IMT
                  WKI (I) = (WKA (I,J,K) - WKA (I -1,J,K))* VIT (I,J,K)*VIT (I -1,J,K)
               END DO
               DO I = 2,IMM
                  TF (I,J,K) = TF (I,J,K) + SOTX (J)* (WKI (I +1) - WKI (I))
!
                  dx(i,j,k,n)=SOTX (J)* (WKI (I +1) - WKI (I))
!
               END DO
            END DO
         END DO
#ifdef SPMD
!      call exch_boundary(dx(1,1,1,n),km,1,1)
#endif
#else
!!$OMP PARALLEL DO PRIVATE (K,J,WKI)
         DO K = 1,KM
            DO J = JSM,JEM
               WKI (1)= 0.0
               DO I = 2,IMT
                  WKI (I) = (ATB (I,J,K,N) - ATB (I -1,J,K,N))* VIT (I,J,K)* VIT (I -1,J,K)
               END DO
               DO I = 2,IMM
                  TF (I,J,K) = TF (I,J,K) + AH3 (I,J,K)* SOTX (J)* (WKI (I +1) - WKI (I))
!
                  dx(i,j,k,n)= AH3 (I,J,K)* SOTX (J)* (WKI (I +1) - WKI (I))
!
               END DO
            END DO
         END DO

#ifdef DEBUG
         call chk_var2d(tf,"df",1)
#endif
 
#ifdef SPMD
!      call exch_boundary(dx(1,1,1,n),km,1,1)
#endif
#endif
 
!-----------------------------------------------------------------------
!     COMPUTE THE EDDY-DIFFUSION TERM :  MERIDIONAL COMPONENT
!-----------------------------------------------------------------------
 
#if (defined BIHAR)
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,IMT
                  WKA (I,J,K) = 0.0D0
               END DO
            END DO
         END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I,wt1,wt2)
         DO K = 1,KM
         DO J = JSM,JEM
         DO I = 2,IMM
            wt1= (ATB (I,J,K,N)- ATB(I,J-1,K,N))*VIT(I,J,K)* VIT (I,J -1,K)
            wt2= (ATB (I,J+1,K,N)- ATB(I,J,K,N))*VIT(I,J,K)* VIT (I,J +1,K)
            wka (I,J,K) = ah3(i,j,k)*(R2D (J)* wt2 - R2C(J)* wt1)
         END DO
         END DO
         END DO

         IF ( NX_proc == 1) then
         DO K=1,KM
            DO J = JSM,JEM
               WKA (1,J,K)= WKA (IMM,J,K)
               WKA (IMT,J,K)= WKA (2,J,K)
            END DO
        END DO
        END IF
#ifdef SPMD
         call exch_boundary(wka(1,1,1),km)
#endif 
!$OMP PARALLEL DO PRIVATE (K,J,I,wt1,wt2)
         DO K = 1,KM
         DO J = JSM,JEM
         DO I = 2,IMM
            wt1= (WKA (I,J,K) - WKA (I,J -1,K))* VIT (I,J,K)* VIT (I,J -1,K)
            wt2= (WKA (I,J+1,K) - WKA (I,J ,K))* VIT (I,J,K)* VIT (I,J +1,K)
            TF (I,J,K) = TF (I,J,K) + (R2D (J)* wt2 - R2C(J)* wt1)
!
            dy(i,j,k,n)=(R2D (J)* wt2 - R2C(J)* wt1)
!          ddy(i,j,k,n)=2.0*ah3(i,j,k)*(R2A(J)*(WKD(I,J+1,K)-WKD(I,J,K))+R2B(J)*(WKD(I,J,K)-WKD(I,J-1,K)))
!
         END DO
         END DO
         END DO
#ifdef SPMD
!      call exch_boundary(dy(1,1,1,n),km,1,1)
!      call exch_boundary(ddy(1,1,1,n),km,1,1)
#endif

#else
!!$OMP PARALLEL DO PRIVATE (K,J,I,wt1,wt2)
         DO K = 1,KM
         DO J = JSM,JEM
         DO I = 2,IMM
            wt1= (ATB (I,J,K,N) - ATB (I,J -1,K,N))* VIT (I,J,K)* VIT (I,J -1,K)
            wt2= (ATB (I,J +1,K,N) - ATB (I,J,K,N))* VIT (I,J,K)* VIT (I,J +1,K)
            TF (I,J,K) = TF (I,J,K) + ah3(i,j,k)* (R2D (J)* wt2 - R2C(J)* wt1)
!
!         dy(i,j,k,n)=ah3(i,j,k)*(R2D (J)* wt2 - R2C(J)* wt1)
!        ddy(i,j,k,n)=2.0*ah3(i,j,k)*(R2A(J)*(AT(I,J+1,K,N)-AT(I,J,K,N))+R2B(J)*(AT(I,J,K,N)-AT(I,J-1,K,N)))
!
         END DO
         END DO
         END DO
#ifdef DEBUG
         call chk_var2d(tf,"df",1)
#endif
#ifdef SPMD
!      call exch_boundary(dy(1,1,1,n),km,1,1)
!      call exch_boundary(ddy(1,1,1,n),km,1,1)
#endif

#endif
#endif
#endif

!-----------------------------------------------------------------------
!     VERTICAL COMPONENT
!-----------------------------------------------------------------------
 
         IF (N == 1)THEN
#if (defined SOLAR)
!     SOLAR SHORTWAVE PENETRATION
        DO K=2,KM-1
!$OMP PARALLEL DO PRIVATE (J,I,wt1,wt2)
        DO J=JSM,JEM
        DO I=2,IMM
            wt1= SWV (I,J)*pen(k-1)*VIT(I,J,K)
            wt2= SWV (I,J)*pen(k)*VIT(I,J,K+1)
            TF (I,J,K)= TF(I,J,K)+(wt1-wt2)*ODZP(K)
!
!lhl0105           penetrate(i,j,k)=(wt1-wt2)*ODZP(k)
            penetrate(i,j,k)= wt2*ODZP(k)
!
        END DO
        END DO
        END DO
!$OMP PARALLEL DO PRIVATE (J,I,wt1,wt2)
        DO J=JSM,JEM
        DO I=2,IMM
           wt1= SWV (I,J)*pen(1)*VIT(I,J,2)
           wt2= SWV (I,J)*pen(km-1)*VIT(I,J,km)
           TF(I,J,1)=TF(I,J,1)-ODZP(1)*wt1
           TF(I,J,km)=TF(I,J,km)+ODZP(km)*wt2
!
!lhl0105            penetrate(i,j, 1)=-wt1*ODZP( 1)
!lhl0105            penetrate(i,j,km)= wt2*ODZP(km)
            penetrate(i,j, 1)= wt1*ODZP( 1)
            penetrate(i,j,km)= wt2*ODZP(km)
!
        END DO
        END DO
#ifdef DEBUG
         call chk_var2d(vit,"vt",1)
         call chk_var2d(tf,"pe",1)
         call chk_var2d(swv,"sw",1)
#endif
#ifdef SPMD
       call exch_boundary(penetrate,km)
#endif
#endif


!==================================
!From here ie. the code about SOLARCHLORO added by linpf
!==================================

#if (defined SOLARCHLORO)
      DO K=2,KM-1
        DO J=JSM,JEM
        DO I=2,IMM
!           if(k.eq.2.and.chloro(i,j).gt.0.1.and.i.eq.10) &
!           write(*,*)'pen_chl_trace',pen_chl(i,j,k-1)/OD0CP,&
!           pen_chl(i,j,k-1),chloro(i,j)
            wt1= SWV(I,J)*pen_chl(I,J,K-1)*VIT(I,J,K)
            wt2= SWV(I,J)*pen_chl(I,J,K)*VIT(I,J,K+1)
            TF (I,J,K)= TF(I,J,K)+(wt1-wt2)*ODZP(K)
!
            penetrate(I,J,K)=(wt1-wt2)*ODZP(K)
!
        END DO
        END DO
        END DO
!$OMP PARALLEL DO PRIVATE (J,I,wt1,wt2)
        DO J=JSM,JEM
        DO I=2,IMM
           wt1= SWV(I,J)*pen_chl(I,J,1)*VIT(I,J,2)
           wt2= SWV(I,J)*pen_chl(I,J,km-1)*VIT(I,J,km)
           TF(I,J,1)=TF(I,J,1)-ODZP(1)*wt1
           TF(I,J,km)=TF(I,J,km)+ODZP(km)*wt2
!
            penetrate(I,J, 1)=-wt1*ODZP(1)
            penetrate(I,J,km)= wt2*ODZP(km)
!
        END DO
        END DO
#ifdef SPMD
       call exch_boundary(penetrate,km)
#endif
#endif

!==================================
!above code about SOLARCHLORO added by linpf
!==================================
         END IF
 
!     EDDY-DIFFUSION
 
        wt1=0
!$OMP PARALLEL DO PRIVATE (K,J,I,wt1,wt2)
        DO K=2,KM-1
           DO J=JSM,JEM
           DO I=2,IMM
              wt1= WKC(I,J,K-1)*(ATB(I,J,K-1,N)-ATB(I,J,K,N))*ODZT(K)*VIT(I,J,K)
              wt2= WKC(I,J,K)*(ATB(I,J,K,N)-ATB(I,J,K+1,N))*ODZT(K+1)*VIT(I,J,K+1)
              TF (I,J,K)= TF(I,J,K)+ODZP(K)*(wt1-wt2)*(1.0D0-AIDIF)
!
              dz(i,j,k,N)=ODZP(K)*(wt1-wt2)*(1.0D0-AIDIF)
!
           END DO
           END DO
        END DO
!$OMP PARALLEL DO PRIVATE (J,I,wt1,wt2)
       DO J = JSM,JEM
       DO I = 2,IMM
           wt1= WKC(I,J,1)*(ATB(I,J,1,N)-ATB(I,J,2,N))*ODZT(2)*VIT(I,J,2)
           wt2= WKC(I,J,km-1)*(ATB(I,J,km-1,N)-ATB(I,J,km,N))*ODZT(km)*VIT(I,J,km)
           TF(I,J,1)=TF(I,J,1)-ODZP(1)*wt1*(1.0D0-AIDIF)
           TF(I,J,km)=TF(I,J,km)+ODZP(km)*wt2*(1.0D0-AIDIF)
!
              dz(i,j, 1,N)=-ODZP( 1)*wt1*(1.0-AIDIF)
              dz(i,j,km,N)= ODZP(km)*wt2*(1.0-AIDIF)
!
       END DO
       END DO
#ifdef DEBUG
         call chk_var2d(tf,"df",1)
#endif
#ifdef SPMD
       call exch_boundary(dz(1,1,1,n),km)
#endif
 
!-----------------------------------------------------------------------
!     SET NEWTONIAN SURFACE BOUNDARY CONDITION
!-----------------------------------------------------------------------
 
         IF (N == 2)THEN
 
!$OMP PARALLEL DO PRIVATE (J,I)
            DO J = JSM,JEM
               DO I = 2,IMM
                  IF (ITNU (I,J) > 0)THEN
#ifdef COUP
!                     STF (I,J) = SSF(I,J)
                     STF (I,J) = SSF(I,J)/ODZP(1)
#else
#ifdef FRC_CORE
                     STF (I,J) = (fresh(i,j)*35.0/1000.0/1000.0&
                   +GAMMA*(SSS(I,J)-ATB (I,J,1,2))*seaice(i,j)/ODZP(1)&
                   +GAMMA*(SSS(I,J)-ATB (I,J,1,2))/ODZP(1)/12.*(1.0-seaice(i,j)))
#else
!                     STF (I,J) = GAMMA * (SSS (I,J) - ATB (I,J,1,2))
                      STF (I,J) = GAMMA * (SSS (I,J) - ATB (I,J,1,2))/ODZP(1)
#endif
#endif
!                     TF (I,J,1) = TF (I,J,1) + STF (I,J)* (1.0- AIDIF)
                     TF (I,J,1) = TF (I,J,1) + STF (I,J)* (1.0- AIDIF)*ODZP(1)
!
!                     NET (I,J,2) = STF (I,J)*(1.0-AIDIF)
                     NET (I,J,2) = STF (I,J)*ODZP(1)
!
                  END IF
               END DO
            END DO
!
        if (nx_proc==1) then
             DO J = JSM,JEM
                      NET (1,J,2) =  NET (IMM,J,2)
                      NET (IMT,J,2) =  NET (2,J,2)
             END DO
        end if
#ifdef SPMD
         call exch_boundary(NET(1,1,2),1,1) 
#endif
        ELSE
!$OMP PARALLEL DO PRIVATE (J,I)
            DO J = JSM,JEM
               DO I = 2,IMM
                IF (ITNU (I,J) > 0)THEN
#ifdef COUP
                 STF (I,J) = TSF(I,J)
#else
#ifdef FRC_CORE
                 STF (I,J) = ((SWV(I,J)+NSWV(I,J))*OD0CP+SEAICE(I,J)*GAMMA*(SST(I,J)-ATB(I,J,1,1))/ODZP(1))
#else
                 STF (I,J) = (SWV(I,J)+NSWV(I,J)-DQDT(I,J)*(SST (I,J) - ATB (I,J,1,1)))*OD0CP
#endif
#endif
                 TF (I,J,1) = TF (I,J,1) + STF (I,J)* ODZP (1)* (1.0D0- AIDIF)
!
!                 NET (I,J,1) = STF (I,J)*ODZP(1)*(1.0-AIDIF)
                 NET (I,J,1) = STF (I,J)*ODZP(1)
!
                END IF
               END DO
            END DO
        if (nx_proc==1) then
             DO J = JSM,JEM
                      NET (1,J,1) =  NET (IMM,J,1)
                      NET (IMT,J,1) =  NET (2,J,1)
             END DO
        end if
#ifdef SPMD
         call exch_boundary(NET(1,1,1),1)
#endif
!       allocate(buffer(imt_global,jmt_global))
!       call local_to_global_4d_double(net(1,1,1),buffer,1,1)
!       if (mytid.eq.0) then
!       write(134,*) buffer
!       close(134)
!       endif
!       call local_to_global_4d_double(swv,buffer,1,1)
!       if (mytid.eq.0) then
!       write(135,*) buffer
!       close(135)
!       endif
!       call local_to_global_4d_double(nswv,buffer,1,1)
!       if (mytid.eq.0) then
!       write(136,*) buffer
!       close(136)
!       endif
!       call local_to_global_4d_double(seaice,buffer,1,1)
!       if (mytid.eq.0) then
!       write(137,*) buffer
!       close(137)
!       endif
!       deallocate(buffer)
!
         END IF 
#ifdef DEBUG
        call chk_var2d(tf,'fi',1)
#endif

      if (nx_proc.ne.1) then
      n2=mod(mytid,nx_proc)
      do j=1,jmt
        if (n2==0) then
         net(1,j,1)=0.0
         net(1,j,2)=0.0
        end if
        if (n2==nx_proc-1)then
         net(imt,j,1)=0.0
         net(imt,j,2)=0.0
        end if
      end do
      end if
 
#if (defined BOUNDARY)
!-----------------------------------------------------------------------
!    boundary condition
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 2,KM
#ifdef SPMD
         do j=1,jmt
            if (j_global(j)>=(jst_global+1).and.j_global(j)<=(jst_global+50)) then
#else
            DO J = JSM,JSM +50
#endif
               DO I = 2,IMM
                  TF (I,J,K)= TF (I,J,K) + VIT (I,J,K)* (RESTORE (I,J,K,&
                             N) - ATB (I,J,K,N))*LAMDA
               END DO
#ifdef SPMD
             endif
#endif
            END DO
         END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 2,KM
#ifdef SPMD
         do j=1,jmt
            if (j_global(j)<=(jmt_global-1).and.j_global(j)>=(jmt_global-50)) then
#else
            DO J = JEM -50,JEM
#endif
               DO I = 2,IMM
                  TF (I,J,K)= TF (I,J,K) + VIT (I,J,K)* (RESTORE (I,J,K,&
                             N) - ATB (I,J,K,N))*LAMDA
               END DO
#ifdef SPMD
             endif
#endif
            END DO
         END DO
 
#endif
!-----------------------------------------------------------------------
!     SOLVE FOR "TAU+1" TRACER AT CENTER OF "T" CELLS
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 2,IMM
!XC
#if (defined TSPAS)
                  VTL (I,J,K) = AT (I,J,K,N) + DTS * TF (I,J,K)
#else
                  VTL (I,J,K) = ATB (I,J,K,N) + C2DTTS * TF (I,J,K)
#endif
!XC
               END DO
            END DO
         END DO
 
!-----------------------------------------------------------------------
!     ADD DT/DT COMPONENT DUE TO IMPLICIT VERTICAL DIFFUSION
!-----------------------------------------------------------------------
 
#if (defined ISO)
         CALL INVTRI (VTL,STF,WKC,AIDIF,C2DTTS)
#endif
 
!---------------------------------------------------------------------
!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!---------------------------------------------------------------------
 
!        if(N.eq.1) then
!	write(*,*) '1'
!	write(*,'(10d15.5)')(vtl(1,j,1),j=1,20)
!	write(*,'(10d15.5)')(vtl(361,j,1),j=1,20)
!        endif
         if (nx_proc==1) then
!$OMP PARALLEL DO PRIVATE (K,J)
         DO K = 1,KM
            DO J = JSM,JEM
               VTL (1,J,K) = VTL (IMM,J,K)
               VTL (IMT,J,K) = VTL (2,J,K)
            END DO
         END DO   
         end if

       call exch_boundary(vtl(1,1,1),km)

!        if(N.eq.1) then
!	write(*,*) '2'
!	write(*,'(10d15.5)')(vtl(1,j,1),j=1,20)
!	write(*,'(10d15.5)')(vtl(361,j,1),j=1,20)
!        endif
               
         if (mod(ist,180) == 1) then
            CALL SMTS (VTL,VIT,KM,fil_lat2)
!           CALL FILTER_TRACER(VTL,VIT_1D,fil_lat2,KM)
         else
             DO J=JSM,JEM
             DO I=1,IMT  
#if (defined TSPAS)
                  VTL (I,J,1) = VTL(I,J,1)-AT (I,J,1,N) - NET(I,J,N)*DTS
#else
                  VTL (I,J,1) = VTL(I,J,1)- ATB (I,J,1,N) - NET(I,J,N)*C2DTTS
#endif
             END DO
             END DO
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
            DO K=2,KM
            DO J=JSM,JEM
            DO I=1,IMT  
#if (defined TSPAS)
                  VTL (I,J,K) = VTL(I,J,K)-AT (I,J,K,N)
#else
                  VTL (I,J,K) = VTL(I,J,K)- ATB (I,J,K,N)
#endif
            END DO
            END DO
            END DO
                  
           CALL SMTS (VTL,VIT,KM,fil_lat1)
!          CALL FILTER_TRACER(VTL,VIT_1D,fil_lat1,KM)
!
!$OMP PARALLEL DO PRIVATE (J,I)
            DO J=JSM,JEM
            DO I=1,IMT  
#if (defined TSPAS)
                  VTL (I,J,1) = AT (I,J,1,N) + VTL(I,J,1) + NET(I,J,N)*DTS
#else
                  VTL (I,J,1) = ATB (I,J,1,N) + VTL(I,J,1) + NET(I,J,N)*C2DTTS
#endif
            END DO
            END DO

!$OMP PARALLEL DO PRIVATE (K,J,I)
            DO K=2,KM
            DO J=JSM,JEM
            DO I=1,IMT  
#if (defined TSPAS)
                  VTL (I,J,K) = AT (I,J,K,N) + VTL(I,J,K)
#else
                  VTL (I,J,K) = ATB (I,J,K,N) + VTL(I,J,K)
#endif
            END DO
            END DO
            END DO
        end if  
!
         if (nx_proc==1) then
!$OMP PARALLEL DO PRIVATE (K,J)
         DO K = 1,KM
            DO J = JSM,JEM
               VTL (1,J,K) = VTL (IMM,J,K)
               VTL (IMT,J,K) = VTL (2,J,K)
            END DO
         END DO
         end if
! 
!        if(N.eq.1) then
!        write(*,*) '4'
!        write(*,'(10d15.5)')(vtl(1,j,1),j=1,20)
!        write(*,'(10d15.5)')(vtl(361,j,1),j=1,20)
!        endif
!-----------------------------------------------------------------------
!     SOLVE FOR "TAU+1" TRACER AT CENTER OF "T" CELLS
!-----------------------------------------------------------------------
!Yu
#ifdef SPMD
      call exch_boundary(vtl,km)
#endif
!Yu
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 1,IMT
                 trend(i,j,k,N)=(VTL(I,J,K)-ATB(I,J,K,N))/C2DTTS*VIT(I,J,K)
               END DO
            END DO
         END DO
 
#ifdef SPMD
         call exch_boundary(trend(1,1,1,N),km)
#endif
 
!XC
#if (defined TSPAS)
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = JST,JET
               DO I = 1,IMT
                  ATB(I,J,K,N)=AT (I,J,K,N)
                  AT (I,J,K,N) = VTL (I,J,K)
               END DO
            END DO
         END DO
      END DO
#else
         IF (IST >= 1)THEN
!$OMP PARALLEL DO PRIVATE (K,J,I)
            DO K = 1,KM
               DO J = JST,JET
                  DO I = 1,IMT
                     ATB (I,J,K,N) = AFT2* AT (I,J,K,N) + AFT1* (ATB (I,&
                                    J,K,N) + VTL (I,J, K))
                  END DO
               END DO
            END DO
         END IF
 
!        if(N.eq.1) then
!        write(*,*) '5'
!        write(*,'(10d15.5)')(vtl(1,j,1),j=1,20)
!        write(*,'(10d15.5)')(vtl(361,j,1),j=1,20)
!        endif
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = JST,JET
               DO I = 1,IMT
                  AT (I,J,K,N) = VTL (I,J,K)
               END DO
            END DO
         END DO
      END DO
#endif
!XC

#else
!$OMP PARALLEL DO PRIVATE (K,J,I)
      do k=1,km
      do j=jst,jet
      do i=1,imt
         atb(i,j,k,1)=12.0D0
         atb(i,j,k,2)=0.0D0
      end do
      end do
      end do
!$OMP PARALLEL DO PRIVATE (N,K,J,I)
      DO N = 1,NTRA
      DO K=1,KM
       DO J=JST,JET
        DO I=1,IMT
         AT(I,J,K,N) = ATB(I,J,K,N)
        ENDDO
       ENDDO
      ENDDO
      END DO
#endif
 
      IST = IST +1
 
#ifdef ISO
      deallocate(K1,K2,K3,adv_vetiso,adv_vbtiso,adv_vntiso)
#endif
      deallocate(stf,tf)
      deallocate(wkb,wkc,wkd)
      deallocate(rit)
      RETURN
      END SUBROUTINE TRACER
 
