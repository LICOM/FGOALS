!  CVS: $Id: work_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module work_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!
!
!     ------------------------------------------------------------------
!     Working Arrays
!     ------------------------------------------------------------------
      real(r8),dimension(imt,jmt):: PXB,PYB,PAX,PAY,WHX,WHY,WGP
      real(r8),dimension(imt,jmt,km)::wka
      real(r8),dimension(:,:,:),allocatable:: work_1,work_2,work_3,temp
!lhl1204      real(r8),dimension(:,:,:),allocatable:: riu
      real(r8),dimension(:,:,:),allocatable:: tmp1,tmp2,uk,vk
      real(r8),dimension(imt,jmt) :: work,work1,work2
      real(r8):: wki(imt),wkj(jmt),wkk(kmp1)
      real(r8),dimension(:,:,:),allocatable:: WKB,WKC,WKD,TF
      real(r8),dimension(:,:),allocatable::stf
      real(r4),dimension(:,:),allocatable::buffer_real4

#ifdef SPMD
      real(r8):: wkj_global(jmt_global)
      real(r8),dimension(:,:), allocatable :: work_global
#endif
!
!
end module work_mod
