!  CVS: $Id: isopyi.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     =================
      SUBROUTINE ISOPYI
!     =================
 
!=======================================================================
!     Initialization for isopycnal mixing scheme
 
!     Redi/Cox version + Gent_McWilliams version
 
!     The mixing unit in MOM2 is cm**2/sec. but in this model, the unit
!     m**2/sec
 
!=======================================================================
use precision_mod 
use param_mod
use pconst_mod
use isopyc_mod
#ifdef SPMD
use msg_mod
#endif
 
      IMPLICIT NONE
      REAL(r8)    :: zdzk,t1,t2,dptmid,xx
 
#ifdef SPMD
      if (mytid==0)then
#endif
      write(6,*)"Beginning-----ISOPYI!"
#ifdef SPMD
      endif 
#endif
!-----------------------------------------------------------------------
!     INPUT SOME PARAMETERS OF IAP MODEL
!-----------------------------------------------------------------------
      zdzk = 0.0D0
 
!!$OMP PARALLEL DO PRIVATE (k),reduction(+:zdzk)
      DO k = 1,km
         zdzk = zdzk + DZP (k)
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k)
      DO k = 1,km
         zt (k)= abs (ZKT (k))
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k)
      DO k = 1,km -1
         dzw (k)= zt (k +1) - zt (k)
      END DO
 
 
      dzw (0)= zt (1)
      dzw (km)= zdzk - zt (km)
 
!$OMP PARALLEL DO PRIVATE (k)
      DO k = 0,km
         dzwr (k)= 1.0D0/ dzw (k)
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k)
      DO k = 1,km
         dzr (k)= ODZP (k)
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (j)
      DO j = 1,jmt
         if (dyt(j) > 0 ) then
            dytr (j)= 1.0D0/ DYR (j)
         else
            dytr (j)= 0.0D0
         end if
      END DO
 
!$OMP PARALLEL DO PRIVATE (j)
      DO j = 1,jmt
         if (dyr(j) > 0 ) then
            dyur (j)= 1.0D0/ DYR (j)
         else
            dyur (j)= 0.0D0
         end if
      END DO
 
!$OMP PARALLEL DO PRIVATE (j)
      DO j = 1,jmt
         if (sint(j) > 0 ) then
            cstrdytr (j)= 1.0D0/ (SINT (j)* DYT (j))
         else
            cstrdytr (j)= 0.0D0
         end if
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (j,k,i)
      DO j = 1,jmt
         DO k = 1,km
            DO i = 1,imt
               tmask (i,k,j)= vit (i,j,k)
            END DO
         END DO
      END DO
 
 
!--------------------------------------------------------------------
 
      slmxr = 100.0D0
 
!IAP  ahisop = 1.d7
      ahisop = 1.d3
 
!     define the isopycnal thickness diffusion coefficient
 
!IAP  athkdf = 1.0d7
      athkdf = 1.0d3
 
!     reference pressure level intervals are defined (see "isopyc.h").
!     "dptlim" must have "nrpl+1" elements
 
!     REMARK: the first and the last elements of "dptlim" must be the
!             depth at the top (0cm) and the maximum bottom depth,
!             respectively. Also, the elements of "dptlim" must be in
!             increasing order.
 
      dptlim (1) = 0.0D0
      dptlim (2) = 1000.0D0
      dptlim (3) = 2000.0D0
      dptlim (4) = 3000.0D0
      dptlim (5) = 4000.0D0
      dptlim (6) = zdzk
 
!-----------------------------------------------------------------------
!     determine the isopycnal reference pressure levels for the "t"
!     grid point levels, using the depths at the "t" grid points as the
!     reference depth (pressure)
!-----------------------------------------------------------------------
 
      DO k = 1,km
         MMM : DO m = 2,nrpl +1
            IF (zt (k) > dptlim (m -1) .AND. zt (k) <= dptlim (m)) THEN
               kisrpl (k) = m -1
               EXIT MMM
            END IF
         END DO MMM

         IF (kisrpl (k) < 1 .OR. kisrpl (k) > nrpl) THEN
            WRITE (6,*) kisrpl (k), k,zt,zkt
            STOP 9100
         END IF
      END DO
 
 
!-----------------------------------------------------------------------
!     the indices used in isopycnal mixing indicating the location of
!     the reference pressure levels in the 20-level table of polynomial
!     expansion variables are computed
 
!     REMARK: because the polynomial expansion coefficients are
!             functions of the reference potential temperature and
!             salinity profiles, at the reference pressure level
!             the corresponding potential temperature and salinity
!             values will be used.
!-----------------------------------------------------------------------
 
      DO m = 1,nrpl
         krplin (m) = 0
      END DO
 
 
      DO m = 2,nrpl +1
         dptmid = 0.5D0* (dptlim (m -1) + dptlim (m))
         IF (dptmid <= zt (1)) THEN
            krplin (m -1) = 1
         ELSE IF (dptmid > zt (km)) THEN
            krplin (m -1) = km
         ELSE IF (dptmid > zt (1) .AND. dptmid <= zt (km)) THEN
            KKK: DO k = 2,km
               IF (zt (k) >= dptmid) THEN
                  t1 = zt (k) - dptmid
                  t2 = dptmid- zt (k -1)
                  IF (t1 > t2) THEN
                     krplin (m -1) = k -1
                  ELSE
                     krplin (m -1) = k
                  END IF
                  EXIT KKK
               END IF
            END DO KKK

         END IF
 
         IF (krplin (m -1) < 1 .OR. krplin (m -1) > km) THEN
            WRITE (6,*) krplin (m -1),m -1
            STOP 9110
         END IF
 
      END DO
 
 
!-----------------------------------------------------------------------
!     the isopycnal diffusion coefficient may be a function of depth. in
!     the default configuration, the isopycnal diffusion coefficient is
!     a constant: "fzisop", which multiplies "ahisop", is set to unity.
!     if "ahisop" varies in the vertical, "fzisop" should contain this
!     variation. the value of "ahisop" should be adjusted accordingly.
!-----------------------------------------------------------------------
 
      DO k = 1,km
         fzisop (k) = 1.0D0
      END DO
!
     do j=1,jmt
        xx=asin(sint(j))*180.0D0/PI
        if (xx > 60.0D0 ) then
            if ( xx > 70.0D0) then
                 F3(j)=0.0D0
            else
                 F3(j)=(70.0D0-xx)*0.1D0
            end if
        else
            F3(j)=1.0D0
        end if
     enddo
!
#ifdef SPMD
      if (mytid==0)then
#endif
      write(6,*)"END-----------ISOPYI!"
#ifdef SPMD
      endif 
#endif
 
      RETURN
      END SUBROUTINE ISOPYI
#else
      SUBROUTINE ISOPYI
      RETURN
      END SUBROUTINE ISOPYI
#endif 
