#define LOGMSG()
!write(mytid+600,'(a,3i4)')"SMUV",__LINE__,k,j
!  CVS: $Id: smuvh.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ========================
      SUBROUTINE SMUV (X,Z,KK,fil_lat)
!     ========================
!     1-D zonal smoother

#include <def-undef.h>
use precision_mod
use param_mod
#ifdef SPMD
use msg_mod
#endif
use pconst_mod, only: sinu,pi,torad
      IMPLICIT NONE
      INTEGER :: JFS1,JFS2,JFN1,JFN2,KK, NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,KK),XS (IMT),Z (IMT,JMT,KM)
      INTEGER     :: NN(JMT), MAX_NN

!lhl      fil_lat=55.D0
!
      MAX_NN = 0
      do j =jst,jmt
         if (sinu(j).le.cos(fil_lat*torad)) then
            NN(j) = int(cos(fil_lat*torad)/sinu(j)*1.2D0)
            if (NN(j) .gt. MAX_NN) MAX_NN = NN(J)
         else 
            NN(j) = 0
         endif
      enddo

!!$OMP PARALLEL DO PRIVATE (K,J,xs)
  
      DO NCY = 1,MAX_NN
         do j =jst,jmt
            if (NN(j) .ge. NCY) then       
               DO K = 1,KK
                  DO I = 1,IMT
                     XS (I)= X (I,J,K)* Z (I,J,K)
                  END DO
                  DO I = 2,IMM
                     X (I,J,K)= (0.5D0*XS(I)+0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,K)
                  END DO
               ENDDO
               if (nx_proc == 1) then
                  DO K = 1,KK
                     X (1,J,K)= X (IMM,J,K)
                     X (IMT,J,K)= X (2,J,K)
                  ENDDO
              end if
            endif
         END DO
         if (nx_proc .ne. 1) then
            call exchange_2D_boundary(x,kk,NN,NCY,0)
            call exchange_2D_boundary(x,kk,NN,NCY,1)
         end if
     END DO

      RETURN
      END SUBROUTINE SMUV



!     ========================
      SUBROUTINE SMTS (X,Z,KK,fil_lat)
!     ========================
!     1-D zonal smoother

#include <def-undef.h>
use precision_mod
use param_mod
#ifdef SPMD
use msg_mod
#endif
use pconst_mod,only: sint,PI,torad
      IMPLICIT NONE

      INTEGER :: JFS1,JFS2,JFN1,JFN2,KK, NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,KK),XS (IMT),Z (IMT,JMT,KM)
      INTEGER     :: NN(JMT), MAX_NN


!     fil_lat=66.D0
!lhl      fil_lat=56.D0
!

      MAX_NN = 0
      do j =jst,jmt
         if (sint(j).le.cos(fil_lat*torad)) then
            NN(j) = int(cos(fil_lat*torad)/sint(j)*1.2D0)
            if (NN(j) .gt. MAX_NN) MAX_NN = NN(J)
         else
            NN(j) = 0
         endif
      enddo

!!$OMP PARALLEL DO PRIVATE (K,J,xs)
      DO NCY = 1,MAX_NN
         do j =jst,jmt
            if (NN(j) .ge. NCY) then
               DO K = 1,KK
                  DO I = 1,IMT
                     XS (I)= X (I,J,K)* Z (I,J,K)
                  END DO
                  DO I = 2,IMM
                     X(I,J,K)=(XS(I)*(1.0D0-0.25D0*Z(I-1,J,K)-0.25D0*Z(I+1,J,K)) &
                               +0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,K)
                  END DO
               ENDDO
               if (nx_proc == 1) then
                  DO K = 1,KK
                     X (1,J,K)= X (IMM,J,K)
                     X (IMT,J,K)= X (2,J,K)
                  ENDDO
               endif
            endif
         enddo
         if (nx_proc .ne. 1) then
            call exchange_2D_boundary(x,kk,NN,NCY,0)
            call exchange_2D_boundary(x,kk,NN,NCY,1)
         end if
     END DO


      RETURN
      END SUBROUTINE SMTS



!     ========================
      SUBROUTINE SMZ0 (X,Z,fil_lat)
!     ========================
!     1-D zonal smoother

#include <def-undef.h>
use precision_mod
      use param_mod
#ifdef SPMD
use msg_mod
#endif
use pconst_mod,only:sint,pi,torad
      IMPLICIT NONE

      INTEGER :: JFS1,JFS2,JFN1,JFN2,NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT),XS (IMT),Z (IMT,JMT,KM)
      INTEGER    :: NN(JMT), MAX_NN


!lhl      fil_lat=54.D0
!
         
      MAX_NN = 0
      do j =jst,jmt
         if (sint(j).le.cos(fil_lat*torad)) then
            NN(j) = int(cos(fil_lat*torad)/sint(j)*1.2D0)
            if (NN(j) .gt. MAX_NN) MAX_NN = NN(J)
         else
            NN(j) = 0
         endif
      enddo

!!$OMP PARALLEL DO PRIVATE (J,nn,xs)
      DO NCY = 1,MAX_NN
         do j =jst,jmt
            if (NN(j) .ge. NCY) then
               DO I = 1,IMT
                  XS (I)= X (I,J)* Z (I,J,1)
               END DO
               DO I = 2,IMM
                  X(I,J)=(XS(I)*(1.0D0-0.25D0*Z(I-1,J,1)-0.25D0*Z(I+1,J,1)) &
                            +0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,1)
               END DO
               if (nx_proc == 1) then
                  X (1,J)= X (IMM,J)
                  X (IMT,J)= X (2,J)
               endif
            endif
         enddo
         if (nx_proc .ne. 1) then
            call exchange_2D_boundary(x,1,NN,NCY,0)
            call exchange_2D_boundary(x,1,NN,NCY,1)
         end if
      END DO


   RETURN
   END SUBROUTINE SMZ0




      SUBROUTINE FILTER_TRACER(X,Z,FIL_LAT,kk)
!
!     1-D zonal smoother

#include <def-undef.h>
use precision_mod
      use param_mod
#ifdef SPMD
use msg_mod
#endif
use pconst_mod,only:sint,pi,torad
      IMPLICIT NONE

      INTEGER :: JFS1,JFS2,JFN1,JFN2,NN, NCY,kk,ix,iy
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,kk),XS (IMT_global,jmt,kk),Z (IMT_global,JMT,kk), tmp(imt_global)
!
      ix = 0
#ifdef SPMD
      if (nx_proc /= 1) then
         call gather_1d(x,xs,kk)
         ix=mod(mytid,nx_proc)
         iy=(mytid-ix)/nx_proc
      end if
#else
      do k=1,kk
      do j=1,jmt
      do i=i,jmt
         xs(i,j,k)=x(i,j,k)
      end do
      end do
      end do
#endif
!
         if ( ix == 0 ) then
         do k=1,kk
         do j =jst , jmt
            if (sint(j).le.cos(fil_lat*torad))then
               nn=int(cos(fil_lat*torad)/sint(j)*1.2D0)
               
               DO NCY = 1,NN
                  DO I = 1,IMT_global
                     tmp (I)= XS(I,J,k)* Z (I,J,k)
                  END DO
               
               DO I = 2,imt_global-1
                  xs(I,J,k)=(tmp(I)*(1.0D0-0.25D0*Z(I-1,J,1)-0.25D0*Z(I+1,J,1)) &
                            +0.25D0*(tmp(I-1)+tmp(I+1)))*Z(I,J,1)
                            
               END DO
                  XS(1,J,K)= XS(imt_global-1,J,K)
                  XS(IMT_global,J,K)= XS(2,J,K)
               END DO
 
            end if
         END DO
         END DO
         end if

#ifdef SPMD
      if (nx_proc /= 1) then
         call distribute_1d(x,xs,kk)
      end if
#endif

      END SUBROUTINE FILTER_TRACER
!
      SUBROUTINE FILTER_MOMENTUM(X,Z,FIL_LAT,kk)
!
!     1-D zonal smoother

#include <def-undef.h>
use precision_mod
      use param_mod
#ifdef SPMD
use msg_mod
#endif
use pconst_mod,only:sinu,pi,torad
      IMPLICIT NONE

      INTEGER :: JFS1,JFS2,JFN1,JFN2,NN, NCY,kk,ix,iy
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,kk),XS (IMT_global,jmt,kk),Z (IMT_global,JMT,kk), tmp(imt_global)
!
      ix = 0 
#ifdef SPMD
      if (nx_proc /= 1) then
         call gather_1d(x,xs,kk)
         ix=mod(mytid,nx_proc)
         iy=(mytid-ix)/nx_proc
      end if
#endif
!
         if (ix == 0 ) then
         do k=1,kk
         do j =jst , jmt
            if (sinu(j).le.cos((fil_lat)*torad))then
               nn=int(cos((fil_lat)*torad)/sinu(j)*1.2D0)
               
               DO NCY = 1,NN
                  DO I = 1,IMT
                     tmp (I)= XS(I,J,k)* Z (I,J,k)
                  END DO
               
               DO I = 2,imt_global-1
                  xs(I,J,K)= (0.5D0*tmp(I)+0.25D0*(tmp(I-1)+tmp(I+1)))*Z(I,J,K)
               END DO
                  XS(1,J,K)= XS(imt_global-1,J,K)
                  XS(imt_global,J,K)= XS(2,J,K)
               END DO
!
            end if
         END DO
         END DO
         end if

#ifdef SPMD
      if (nx_proc /= 1) then
         call distribute_1d(x,xs,kk)
      end if
#endif

      END SUBROUTINE FILTER_MOMENTUM
