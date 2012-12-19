!  CVS: $Id: smag.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!      this program is Smagrinsky horizontal viscosity
!      based on MOM2 Manual and Rosati A. & K. Miyakoda (1988)
!      writted by liu hailong jun 2001.
 
#include <def-undef.h>
#if ( defined SMAG)
      SUBROUTINE smag2 (KK)
use precision_mod 
use param_mod
use pconst_mod
use work_mod
 
      IMPLICIT none
      INTEGER :: KK
      REAL(r8)    :: abcd,bond
 
#if (defined SMAG_FZ)
!$OMP PARALLEL DO PRIVATE (J,I)
      DO j = 2,jmm
         DO i = 2,imm
            wka (i,j,13)= ((0.5D0*oux(j)*(-wka(i-1,j,11)+wka(i+1,j,11)) &
                          - (r1e (j)* wka (i,j,12) - r1f (j)* wka (i,j-1,12) &
                          +r1e(j+1)*wka(i,j+1,12)-r1f(j+1)*wka(i,j,12))))
            wka (i,j,14)= -((0.5D0*oux(j)*(-wka(i-1,j,12)+wka(i+1,j,12)) &
                          +(r1e(j)*wka(i,j,11)-r1f(j)* wka(i,j-1,11)  &
                          + r1e(j+1)*wka (i,j+1,11)-r1f(j+1)*wka(i,j,11))))
         END DO
      END DO
!
     if (nx_proc == 1) then 
     do j =2,jmm
         wka (1,j,13)= wka (imm,j,13)
         wka (1,j,14)= wka (imm,j,14)
         wka (imt,j,13)= wka (2,j,13)
         wka (imt,j,14)= wka (2,j,14)
     end do
     end if
 
!$OMP PARALLEL DO PRIVATE (I)
      DO i = 1,imt
         wka (i,1,13)= 0.0D0
         wka (i,1,14)= 0.0D0
         wka (i,jmt,13)= 0.0D0
         wka (i,jmt,14)= 0.0D0
      END DO
 
!-new
#else
!$OMP PARALLEL DO PRIVATE (I)
      DO i = 1,imt
         wka (i,1,13)= 0.0D0
         wka (i,1,14)= 0.0D0
         wka (i,jmt,13)= 0.0D0
         wka (i,jmt,14)= 0.0D0
      END DO
!
!$OMP PARALLEL DO PRIVATE (J,I)
      DO j = 2,jmm
         DO i = 2,imm
            wka (i,j,13)= (0.5D0* oux (j)* (wka (i +1,j,11) - wka (i -1,j,11)) &
                        - r1e(j)*wka(i,j +1,12) + r1f (j)* wka (i,j -1,12))
            wka (i,j,14)= - (0.5D0* oux (j)*(wka(i+1,j,12)-wka(i-1,j,12)) &
                          + r1e (j)* wka (i,j +1,11) - r1f (j)* wka (i,j -1,11))
         END DO
      END DO
!
     if (nx_proc == 1) then
         do j=2,jmm
         wka (1,j,13)= wka (imm,j,13)
         wka (1,j,14)= wka (imm,j,14)
         wka (imt,j,13)= wka (2,j,13)
         wka (imt,j,14)= wka (2,j,14)
         end do
    end if
 
 
#endif
 
!
#ifdef  SPMD 
       call exch_boundary(wka(1,1,13),1)
       call exch_boundary(wka(1,1,14),1)
#endif
!
!$OMP PARALLEL DO PRIVATE (J,I)
      DO j = jst,jmt
         DO i = 1,imt
!lhl            wka (i,j,15)= sqrt (wka (i,j,13)**2+ wka (i,j,14)**2)
            wka (i,j,15)= sqrt (2D0*wka (i,j,13)**2+ 2D0*wka (i,j,14)**2)
         END DO
      END DO
 
 
!lhl!!$OMP PARALLEL DO PRIVATE (j,i)
!lhl      DO j = jst,jmt
!lhl         DO i = 1,imt
!lhl            am3 (i,j,KK)= am (j)
!lhl         END DO
!lhl      END DO
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO j = 1,jmt
         DO i = 1,imt
            amx (i,j,KK)=  wka (i,j,15)* cxu (j)       
            amy (i,j,KK)=  wka (i,j,15)* cyu (j)       
!lhl            am3 (i,j,KK)= am3 (i,j,kk) + wka (i,j,15)* cxu (j)       
         END DO
      END DO

!#ifdef  SPMD 
!	if (mytid == 0) then
!	write(111,*)amx
!        endif
!#endif
 
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO j = jst,jmt
         DO i = 1,imt
!lhl            wka (i,j, 7)= am3 (i,j,KK)* wka (i,j,13)* viv (i,j,KK)
!lhl            wka (i,j, 8)= am3 (i,j,KK)* wka (i,j,14)* viv (i,j,KK)
!lhl            wka (i,j, 9)= am3 (i,j,KK)* wka (i,j,14)* viv (i,j,KK)
!lhl            wka (i,j,10)= - am3 (i,j,KK)* wka (i,j,13)* viv (i,j,KK)
            wka (i,j, 7)= amx (i,j,KK)* wka (i,j,13)* viv (i,j,KK)
            wka (i,j, 8)= amy (i,j,KK)* wka (i,j,14)* viv (i,j,KK)
            wka (i,j, 9)= amx (i,j,KK)* wka (i,j,14)* viv (i,j,KK)
            wka (i,j,10)= - amy (i,j,KK)* wka (i,j,13)* viv (i,j,KK)
         END DO
      END DO
 
 
      END SUBROUTINE smag2
 
 
 
      SUBROUTINE smag3
use precision_mod 
use param_mod
use pconst_mod
use work_mod
use dyn_mod
      IMPLICIT none
 
      REAL(r8)    :: dst (imt,jmt,2),dd (imt,jmt),abcd
 
      DO k = 1,km
#if (defined SMAG_FZ)
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = 2,jmt
            DO i = 2,imm
               dst(i,j,1)=0.5D0*otx(j)*((utl(i,j,k)-utl(i-1,j,k))+(utl(i+1,j,k)-utl(i,j,k))) &
                         -(r1e(j)*vtl(i,j,k)-r1f(j)*vtl(i,j-1,k)+r1e(j)*vtl(i,j+1,k)       &
                         -r1f(j+1)*vtl(i,j,k))
               dst(i,j,2)=-(0.5D0*otx(j)*((vtl(i,j,k)-vtl(i-1,j,k))+(vtl(i+1,j,k)-vtl(i,j,k))) &
                          +(r1e(j)*utl(i,j,k)-r1f(j)*utl(i,j-1,k)+r1e(j+1)*utl(i,j+1,k)      &
                          -r1f(j+1)*utl(i,j,k)))
            END DO
         END DO
      if (nx_proc==1) then
         do j=2,jmt
            dst (1,j,1)= dst (imm,j,1)
            dst (1,j,2)= dst (imm,j,2)
            dst (imt,j,1)= dst (2,j,1)
            dst (imt,j,2)= dst (2,j,2)
         end do
     end if
 
         DO i = 1,imt
            dst (i,1,1)= 0.0D0
            dst (i,1,2)= 0.0D0
            dst (i,jmt,1)= 0.0D0
            dst (i,jmt,2)= 0.0D0
         END DO
#ifdef  SPMD                                                   
       call exch_boundary(dst(1,1,1),1)
       call exch_boundary(dst(1,1,2),1)
#endif 
 
#else
         DO i = 1,imt
            dst (i,JMT,1)= dst (i,JMM,1)
            dst (i,JMT,2)= dst (i,JMM,2)
            dst (i,1,1)= dst (i,2,1)
            dst (i,1,2)= dst (i,2,2)
         END DO
!
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = 2,jmm
            DO i = 2,imm
               dst (i,j,1)= viv(i,j,k)*(0.5D0*oux(j)*(utl(i+1,j,k)-utl(i-1,j,k)) &
                           -r1e(j)*vtl(i,j+1,k)+r1f(j)*vtl(i,j-1,k))     
               dst (i,j,2)= - viv(i,j,k)*(0.5D0*oux(j)*(vtl(i+1,j,k)-vtl(i-1,j,k)) &
                            +r1e (j)* utl (i,j +1,k) - r1f (j)* utl (i,j -1,k))   
            END DO
         END DO
!
      if (nx_proc==1) then
         do j=2,jmm
            dst (1,j,1)= dst (2,j,1)
            dst (1,j,2)= dst (2,j,2)
            dst (imt,j,1)= dst (imm,j,1)
            dst (imt,j,2)= dst (imm,j,2)
         end do
      end if 
 
#ifdef  SPMD                                                   
       call exch_boundary(dst(1,1,1),1)
       call exch_boundary(dst(1,1,2),1)
#endif 
 
#endif
 
         DO j = jsm,jmm
            DO i = 1,imt
               dd (i,j)= sqrt (dst (i,j,1)**2+ dst (i,j,2)**2)
            END DO
         END DO
 
 
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst,jmt
            DO i = 1,imt
               am3 (i,j,k)= am (j)
            END DO
         END DO
 
 
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jsm,jmm
            DO i = 1,imt
               am3 (i,j,k)= am3 (i,j,k)*rr+dd (i,j)* cxu(j)*rr    
            END DO
         END DO
!
      END DO
 
     call exch_boundary(am3,km)
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO k = 1,km
         DO j = jsm,jmm
            DO i = 2,imm
               IF (vit (i,j,k) > 0.5) THEN
                  ah3 (i,j,k)= vit (i,j,k)* (am3 (i,j,k)* viv (i,j,k) &
                             + am3 (i +1,j,k)* viv (i +1,j,k) &
                             + am3 (i,j -1,k)* viv (i,j -1,k) &
                             + am3 (i +1,j -1,k)* viv (i +1,j -1,k))/ &
                             (viv (i,j,k) + viv (i +1,j,k) &
                             + viv (i,j -1,k) + viv (i +1,j -1,k))          
               END IF
            END DO
         END DO
      END DO
!
     if (nx_proc==1) then
!$OMP PARALLEL DO PRIVATE (K,J)
         do k=1,km
         do j=jsm,jem
            ah3(1,j,k)=ah3(imm,j,k)
            ah3(imt,j,k)=ah3(2,j,k)
         end do
         end do
     end if
#ifdef SPMD
      call exch_boundary(ah3,km)
#endif 
 
      END SUBROUTINE smag3
 
 
#else
      SUBROUTINE smag2 ()
      END SUBROUTINE smag2
 
 
      SUBROUTINE smag3 ()
      END SUBROUTINE smag3
 
#endif 
