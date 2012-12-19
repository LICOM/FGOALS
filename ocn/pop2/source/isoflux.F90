!  CVS: $Id: isoflux.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     ==============================
      SUBROUTINE ISOFLUX (MTRACE)
!     ==============================
 
!     isopycnal diffusive tracer fluxes are computed.
use precision_mod 
use param_mod
use pconst_mod
use tracer_mod
use isopyc_mod
use work_mod
#ifdef SPMD
use msg_mod
#endif
 
      IMPLICIT NONE
 
      INTEGER :: mtrace
      REAL(r8):: p5,p25,c1,c0,fxa,fxb,fxc,fxe
 
      allocate (work_1(imt,jmt,km),work_2(imt,jmt,km))
      allocate (temp(imt,jmt,km))
      allocate (work_3(imt,jmt,0:km))
!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------
      p5 = 0.5D0
 
      c0 = 0.0D0
      c1 = 1.0D0
      p25 = 0.25D0
 
      m = mtrace
 
!-----------------------------------------------------------------------
!     first compute the vertical tracer flux "temp" at the northern
!     face of "t" cells.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 2,km -1
         DO j = 2,jmm
            DO i = 1,imt
               temp (i,j,k)= p25* dzr (k)* (atb (i,j +1,k -1,m) - atb ( &
                  i,j +1,k +1,m) &
                   + atb (i,j,k -1,m) - atb (i,j, k +1,m))    
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     now consider the top level, assuming that the surface tracer
!     values are the same as the ones at "k"=1
!-----------------------------------------------------------------------
 
      k = 1
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 2,jmm
         DO i = 1,imt
            temp (i,j,k) = 0.25* dzr (k)* (atb (i,j +1,k,m) - atb (i,   &
                          j +1,k +1,m)+atb (i,j,k,m) - atb (i,j, k +1,m))
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     finally, consider the bottom level. the extrapolative estimator
!     is used to compute the tracer values at the ocean bottom.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 1,jmt
         DO i = 1,imt
            temp (i,j,km)= 0.0D0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (j,i,k,fxa,fxb,fxc,fxe)
      DO j = 2,jmm
         DO i = 1,imt
            k = min (ITNU (i,j),ITNU (i,j +1))
            IF (k /= 0) THEN
               fxe = dzw (k -1) + dzw (k)
               fxa = 0.5D0* (atb (i,j +1,k -1,m) + atb (i,j,k -1,m))
               fxb = 0.5D0* (atb (i,j +1,k,m) + atb (i,j,k,m))
               fxc = dzwr (k -1)* (fxb * fxe- fxa * dzw (k))
               temp (i,j,k) = dzr (k)* (0.5D0* (fxa + fxb) - fxc)
            END IF
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     compute of meridional tracer flux
!     first calculate the effects of purely horizontal diffusion along
!     isopycnal, the background horizontal diffusion has been computed
!     before called this subroutine. (jxz)
!     add in the effects of the along isopycnal diffusion computed
!     using "K2" component of the tensor and apply land/sea masks
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO k = 1,km
         DO i = 1,imt
            work_1 (i,1,k)= 0.0D0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 1,imt
               work_1 (i,j,k)= ( ahisop * dyur (j)* (atb (i,j +1,k,m)- atb (i,j,k,m)) + &
               ahisop * K2 (i,k,j,3)* temp (i,j,k) )* &
               SINU (j)* vit (i,j,k)* vit (i,j +1,k) 
!
               ddy_iso(i,j,k,m)=work_1(i,j,k)/SINU(J)
!
            END DO
         END DO
      END DO
#ifdef SPMD
      call exch_boundary(work_1,km)
      call exch_boundary(ddy_iso(1,1,1,m),km)
#endif
 
 
!-----------------------------------------------------------------------
!     compute the vertical tracer flux "temp" at the eastern
!     face of "t" cells.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 2,km -1
         DO j = 2,jmm
            DO i = 1,imm
               temp (i,j,k)= p25* dzr (k)* (atb (i +1,j,k -1,m) - atb ( &
                             i+1,j,k+1,m)+atb(i,j,k-1,m)-atb(i,j,k+1,m))     
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     now consider the top level, assuming that the surface tracer
!     values are the same as the ones at "k"=1
!-----------------------------------------------------------------------
 
      k = 1
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 2,jmm
         DO i = 1,imm
            temp (i,j,k)= p25* dzr (k)* (atb (i +1,j,k,m) - atb (i +1,j,&
               k +1,m) + atb (i,j,k,m) - atb (i,j,k +1,m))  
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     finally, consider the bottom level. the extrapolative estimator
!     is used to compute the tracer values at the ocean bottom.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 2,jmm
         DO i = 1,imm
            temp (i,j,km) = 0.0D0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (j,i,k,fxa,fxb,fxc,fxe)
      DO j = 2,jmm
         DO i = 1,imm
            k = min (ITNU (i,j),ITNU (i +1,j))
            IF (k /= 0) THEN
               fxe = dzw (k -1) + dzw (k)
               fxa = p5* (atb (i,j,k -1,m) + atb (i +1,j,k -1,m))
               fxb = p5* (atb (i,j,k,m) + atb (i +1,j,k,m))
               fxc = dzwr (k -1)* (fxb * fxe- fxa * dzw (k))
               temp (i,j,k) = dzr (k)* (p5* (fxa + fxb) - fxc)
            END IF
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     compute of zonal tracer flux
!     first calculate the effects of purely horizontal diffusion along
!     isopycnal, the background horizontal diffusion has been computed
!     before called this subroutine. (jxz)
!     add in the effects of the along isopycnal diffusion computed
!     using "K1" component of the tensor and apply land/sea masks
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 1,imm
               work_2 (i,j,k)= ( ahisop * OTX (j)* (atb (i +1,j,k,m)    &
                              - atb (i,j,k,m)) + &
               ahisop*K1(i,k,j,3)*temp(i,j,k) )*vit(i+1,j,k)* vit(i,j,k) 
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     compute the vertical tracer flux "work_3" containing the K31
!     and K32 components which are to be solved explicitly. The K33
!     component will be treated semi-implicitly
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 2,km
         DO j = 2,jmm
            DO i = 2,imm
               work_3 (i,j,k -1) = ahisop * p25* vit (i,j,k)* ( &
               OTX (j)* K3 (i,k -1,j,1)* &
               (vit (i -1,j,k )* (atb (i,j,k,m) - atb (i -1,j,k,m)) &
               + vit (i -1,j,k -1)* (atb (i,j,k -1,m) - atb (i -1,j,k -1,m)) &
               + vit (i +1,j,k )* (atb (i +1,j,k,m) - atb (i,j,k,m)) &
               + vit (i +1,j,k -1)* (atb (i +1,j,k -1,m) - atb (i,j,    &
                                  k -1,m))) + &
               dytr (j)* K3 (i,k -1,j,2)* &
               (vit (i,j -1,k )* (atb (i,j,k,m) - atb (i,j -1,k,m)) &
               + vit (i,j -1,k -1)* (atb (i,j,k -1,m) - atb (i,j -1,k -1,m)) &
               + vit (i,j +1,k )* (atb (i,j +1,k,m) - atb (i,j,k,m)) &
                                   + vit (i,j +1,k -1)* (atb (i,j +1,   &
                                    k -1,m) - atb (i,j,k -1,m))) )  
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     at ocean surface the flux is set to zero to reflect the no tracer
!     flux condition. Same condition is also imposed at ocean bottom.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 2,jmm
         DO i = 2,imm
            work_3 (i,j, 0)= 0.0D0
            work_3 (i,j,km)= 0.0D0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 2,imm
               tf (i,j,k) = tf (i,j,k) &
               + cstrdytr (j)* (work_1 (i,j,k) - work_1 (i,j -1,k)) &
               + OTX (j) * (work_2 (i,j,k) - work_2 (i -1,j,k)) &
               + dzr (k) * (work_3 (i,j,k -1) - work_3 (i,j,k))
!
               dx_iso(i,j,k,m)= OTX (j) * (work_2 (i,j,k) - work_2 (i -1,j,k)) 
               dy_iso(i,j,k,m)= cstrdytr (j)* (work_1 (i,j,k) - work_1 (i,j -1,k))
               dz_iso(i,j,k,m)= dzr (k) * (work_3 (i,j,k -1) - work_3 (i,j,k))
!
            END DO
         END DO
      END DO
#ifdef SPMD
       call exch_boundary(dx_iso(1,1,1,m),km,1,1)
       call exch_boundary(dy_iso(1,1,1,m),km,1,1)
       call exch_boundary(dz_iso(1,1,1,m),km,1,1)
#endif
 
 
!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal velocity mixing
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 2,imm
               work_1 (i,j,k) = adv_vntiso (i,k,j)* (atb (i,j +1,k,m)   &
                                + atb (i,j,k,m))
!
               aay_iso(i,j,k,m) =work_1(i,j,k)*p5/SINU(J)               
!
            END DO
         END DO
      END DO
#ifdef SPMD 
      call exch_boundary(aay_iso(1,1,1,m),km,1,1)
#endif
 
#ifdef SPMD 
!PY      if (mytid==0) then
!PY!$OMP PARALLEL DO PRIVATE (k,i)
!PY      DO k = 1,km
!PY         DO i = 2,imm
!PY            work_1 (i,1,k) = 0.0D0
!PY         END DO
!PY      END DO
!PY      end if

      call exch_boundary(work_1,km,0,1)
#else
!$OMP PARALLEL DO PRIVATE (k,i)
      DO k = 1,km
         DO i = 2,imm
            work_1 (i,1,k) = 0.0D0
         END DO
      END DO
#endif
 
 
!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal velocity mixing
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 1,imm
               work_2 (i,j,k) = adv_vetiso (i,k,j)* (atb (i +1,j,k,m)   &
                                + atb (i,j,k,m))
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     compute the vertical component of the isopycnal velocity mixing
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 2,jmm
         DO i = 2,imm
            work_3 (i,j, 0)= 0.0D0
            work_3 (i,j,km)= 0.0D0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 2,km
         DO j = 2,jmm
            DO i = 2,imm
               work_3 (i,j,k -1)= adv_vbtiso (i,k -1,j)* (atb (i,j,k,m) &
                                + atb (i,j,k -1,m))
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 2,jmm
            DO i = 2,imm
               tf (i,j,k) = tf (i,j,k) &
               - p5* cstrdytr (j)* (work_1 (i,j,k) - work_1 (i,j -1,k)) &
               - p5* OTX (j) * (work_2 (i,j,k) - work_2 (i -1,j,k)) &
               - p5* dzr (k) * (work_3 (i,j,k -1) - work_3 (i,j,k)) 
!
              ay_iso(i,j,k,m)=- p5* cstrdytr (j)* (work_1 (i,j,k) - work_1 (i,j -1,k)) 
              ax_iso(i,j,k,m)=- p5* OTX (j) * (work_2 (i,j,k) - work_2 (i -1,j,k)) 
              az_iso(i,j,k,m)=- p5* dzr (k) * (work_3 (i,j,k -1) - work_3 (i,j,k)) 
!
            END DO
         END DO
      END DO
#ifdef SPMD
       call exch_boundary(ax_iso(1,1,1,m),km)
       call exch_boundary(ay_iso(1,1,1,m),km)
       call exch_boundary(az_iso(1,1,1,m),km)
#endif
      deallocate (work_1,work_2,work_3,temp)
 
      RETURN
      END SUBROuTINE ISOFLUX
 
 
#else
      SUBROUTINE ISOFLUX ()
      RETURN
      END SUBROUTINE ISOFLUX
#endif 
