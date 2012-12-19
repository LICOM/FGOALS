!  CVS: $Id: accumm.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE ACCUMM
!     =================
 
#include <def-undef.h>
use param_mod
use dyn_mod
use tracer_mod
use output_mod
use pconst_mod
use forc_mod, only: su,sv,lthf,sshf,lwv,swv,fresh,runoff
      IMPLICIT NONE

#ifdef COUP
       call exchange_2d(runoff,1,1) !LPF 20120904
       call exchange_2d(fresh,1,1)
       call exchange_2d(lthf,1,1)
       call exchange_2d(sshf,1,1)
       call exchange_2d(lwv,1,1)
#endif

 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = 1,JMT
         DO I = 1,IMT
            Z0MON (I,J)= Z0MON (I,J) + H0 (I,J)
            HIMON (I,J)= HIMON (I,J) + runoff(i,j) !HI (I,J) !LPF 20120824
            HDMON (I,J)= HDMON (I,J) + ALEAD (I,J)* HI (I,J)
!
!linpf091126
            sumon (I,J)= sumon (I,J) + su(I,J)        !U windstress
            svmon (I,J)= svmon (I,J) + sv(I,J)        !V windstress
            lthfmon (I,J)= lthfmon (I,J) + lthf (I,J) !latent flux
            sshfmon (I,J)= sshfmon (I,J) + sshf (I,J) !sensible flux
            lwvmon (I,J)= lwvmon (I,J) + lwv (I,J)    !long wave flux
            swvmon (I,J)= swvmon (I,J) + swv (I,J)    !shortwave flux
!linpf091126 
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               WSMON (I,J,K)= WSMON (I,J,K) + WS (I,J,K)
               TSMON (I,J,K)= TSMON (I,J,K) + AT (I,J,K,1)
               SSMON (I,J,K)= SSMON (I,J,K) + (AT (I,J,K,2)*1000.+35.)
               USMON (I,J,K)= USMON (I,J,K) + U (I,J,K)
               VSMON (I,J,K)= VSMON (I,J,K) - V (I,J,K)
#if (defined SMAG_OUT)
               AM3MON (I,J,K)= AM3MON (I,J,K) + AM3 (I,J,K)
#endif
            END DO
         END DO
      END DO

!$OMP PARALLEL DO PRIVATE (N,K,J,I)
      DO N = 1,NTRA
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               trendmon (I,J,K,N)= trendmon (I,J,K,N)+trend(i,j,k,n)
               axmon (I,J,K,N)= axmon (I,J,K,N)+ax(i,j,k,n)
               aymon (I,J,K,N)= aymon (I,J,K,N)+ay(i,j,k,n)
               azmon (I,J,K,N)= azmon (I,J,K,N)+az(i,j,k,n)
               dxmon (I,J,K,N)= dxmon (I,J,K,N)+dx(i,j,k,n)
               dymon (I,J,K,N)= dymon (I,J,K,N)+dy(i,j,k,n)
               dzmon (I,J,K,N)= dzmon (I,J,K,N)+dz(i,j,k,n)
!
               ddymon (I,J,K,N)= ddymon (I,J,K,N)+ddy(i,j,k,n)
#ifdef ISO
               axmon_iso (I,J,K,N)= axmon_iso (I,J,K,N)+ax_iso(i,j,k,n)
               aymon_iso (I,J,K,N)= aymon_iso (I,J,K,N)+ay_iso(i,j,k,n)
               azmon_iso (I,J,K,N)= azmon_iso (I,J,K,N)+az_iso(i,j,k,n)
               dxmon_iso (I,J,K,N)= dxmon_iso (I,J,K,N)+dx_iso(i,j,k,n)
               dymon_iso (I,J,K,N)= dymon_iso (I,J,K,N)+dy_iso(i,j,k,n)
               dzmon_iso (I,J,K,N)= dzmon_iso (I,J,K,N)+dz_iso(i,j,k,n)
!
               aaymon_iso (I,J,K,N)= aaymon_iso (I,J,K,N)+aay_iso(i,j,k,n)
               ddymon_iso (I,J,K,N)= ddymon_iso (I,J,K,N)+ddy_iso(i,j,k,n)
#endif
            END DO
         END DO
      END DO
      END DO
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               penmon (I,J,K)= penmon (I,J,K)+penetrate(i,j,k)
            END DO
         END DO
      END DO
!
!$OMP PARALLEL DO PRIVATE (N,J,I)
      DO N = 1,NTRA
         DO J = 1,JMT
            DO I = 1,IMT
               netmon (I,J,N)= netmon (I,J,N) + net (I,J,N)
            END DO
         END DO
      END DO
!
!lhl1204
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = 1,JMT
            DO I = 1,IMT
               mldmon (I,J)= mldmon (I,J) + amld (I,J)/100.
            END DO
         END DO
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               akmmon (I,J,K)= akmmon (I,J,K) + akmu(I,J,K)
               aktmon (I,J,K)= aktmon (I,J,K) + akt(I,J,K,1)
               aksmon (I,J,K)= aksmon (I,J,K) + akt(I,J,K,2)
            END DO
         END DO
      END DO
!lhl1204
!
      RETURN
      END SUBROUTINE ACCUMM
 
 
