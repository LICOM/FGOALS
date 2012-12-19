
subroutine cldinti ()

   use shr_kind_mod, only: r8 => shr_kind_r8
   use spmd_utils,   only: masterproc
   use pmgrid,       only: plev, plevp
   use cldconst,     only: k700
   use abortutils,   only: endrun
   use hycoef,       only: hypm 
   use cam_logfile,  only: iulog
   implicit none

   integer :: k

!
! Find vertical level nearest 700 mb
!
   k700 = 1
   do k=1,plev-1
      if (hypm(k) < 7.e4_r8 .and. hypm(k+1) >= 7.e4_r8) then
         if (7.e4_r8-hypm(k) < hypm(k+1)-7.e4_r8) then
            k700 = k
         else
            k700 = k + 1
         end if
         goto 20
      end if
   end do

   call endrun ('CLDINTI: model levels bracketing 700 mb not found')

20 continue

   if (masterproc) then
      write(iulog,*)'CLDINTI: model level nearest 700 mb is',k700,'which is',hypm(k700),'pascals'
   end if

   return
end subroutine cldinti
