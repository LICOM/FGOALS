module sw_mod
#include <def-undef.h> 
!     
use precision_mod
use param_mod  
#if(defined SOLARCHLORO)
      integer,parameter :: M_chl=31
      integer,parameter :: nsub = 400
      INTEGER :: mc,chlindx 
      real(r8) :: arg_max,arg,logchl
      real(r8) :: percm,chlamnt
      real(r8) :: c0,c1,A_1,B_1,A_2,B_2,w11,w12 
      real(r8),dimension(31) :: chloc,A_11,A_12,B_11,B_12
      real(r8) :: chl_temp,chlmin,chlmax,dlogchl 
      real(r8) :: ztr(km),Tr(0:km-1,0:nsub)
#endif
end module sw_mod
