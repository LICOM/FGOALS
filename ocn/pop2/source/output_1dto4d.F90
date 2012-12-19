!  CVS: $Id: output1d.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ============================================
      subroutine output1d()
!     ============================================
! output 1D data by LPF 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use output_mod
use dyn_mod
use tracer_mod
use work_mod
      IMPLICIT NONE
#include <netcdf.inc>
      integer :: itag
!      logical :: hist_output,rest_output 
!

          end subroutine output1d 


!     =================================
     subroutine output2d(start3,count3,ixx,jyy,kk,ncid,iret,x2d_id,x2d,t2_cdf)
!     =================================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
      IMPLICIT NONE
#include <netcdf.inc>
      integer :: kk,x2d_id,ncid,iret,ixx,jyy
      real(r4)    :: t2_cdf(ixx,jyy,kk),x2d(ixx,jyy)
      integer:: start3(3),count3(3)

!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = 1,jyy
            DO i = 1,ixx
                  t2_cdf (i,j,1)= x2d (i,j)/ (nmonth (mon0))
            END DO
         END DO

         iret = nf_put_vara_real (ncid,x2d_id,start3, count3, t2_cdf)
         CALL check_err (iret)

             return
        end subroutine output2d 

