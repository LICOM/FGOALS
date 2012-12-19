!  CVS: $Id: forc_mod.F90,v 1.7 2003/08/25 07:47:52 lhl Exp $
module forc_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!
!
!     ------------------------------------------------------------------
!     Forcing Fields
!     ------------------------------------------------------------------
      real(r8),allocatable,dimension(:,:,:)::su3,sv3,psa3,tsa3,qar3,uva3, &
                   swv3,cld3,sss3,sst3&
!lhl
                    ,nswv3,dqdt3,chloro3&
!lhl
!lhl090730
                   ,wspd3,wspdu3,wspdv3,lwv3,seaice3,rain3,snow3,runoff3
!lhl090730
      real(r8),dimension(imt,jmt)::su,sv,psa,tsa,sss,swv,uva,qar,cld,&
                    ddd,qqq,sst&
!lhl
                   ,nswv,dqdt,chloro&
!lhl
!lhl090730
                   ,lwv,seaice,rain,snow,fresh,runoff&
!lhl090730
!linpf091126
                   ,lthf,sshf !only for output
!linpf091126

!lhl1204
      real(r8),dimension(imt,jmt)::USTAR,BUOYTUR,BUOYSOL
!lhl1204
!
!#ifdef SPMD
      real(r8),allocatable,dimension(:,:,:)::su3_io,sv3_io,psa3_io,tsa3_io,qar3_io,uva3_io, &
                    swv3_io,cld3_io,sss3_io,sst3_io&
                   ,nswv3_io,dqdt3_io,chloro3_io& 
!lhl
!lhl090730
                   ,wspdu3_io,wspdv3_io,lwv3_io,seaice3_io,rain3_io,snow3_io,runoff3_io
!lhl090730
                                                        
#if (defined FRC_DAILY)
      real(r4),dimension(imt_global,jmt_global)::su_in_io,sv_in_io
      real(r8),dimension(imt_global,jmt_global)::su_io,sv_io
#endif
!lhl
!#endif
!
!
#if (defined BOUNDARY)
      real(r8),dimension(imt,jmt,km,ntra)::restore
      ! lihuimin 2012.6.18
      real(r8),dimension(imt,jmt_global,km,ntra)::restore_io
      ! modi end
#endif
!
#ifdef SPMD
       real(r4),dimension(imt_global,jmt_global):: tsf_global,ssf_global,su_global,sv_global,swv_global
#endif

       real(r8),dimension(imt,jmt):: tsf,ssf
end module forc_mod
