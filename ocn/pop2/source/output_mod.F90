!  CVS: $Id: output_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module output_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     Output Arrays
!     ------------------------------------------------------------------
      real(r4),dimension(imt,jmt)::z0mon,himon,hdmon
!
!linpf091126      
      real(r4),dimension(imt,jmt)::lthfmon,sshfmon,lwvmon,swvmon
      real(r4),dimension(imt,jmt)::sumon,svmon
!linpf091126
!
      real(r4),dimension(imt,jmt,km)::wsmon,tsmon,ssmon,usmon,vsmon
      real(r4),dimension(imt,jmt,2):: icmon
      real(r4),dimension(imt,jmt,NTRA)::netmon
!
!lhl1204
      real(r4),dimension(imt,jmt)::mldmon
      real(r4),dimension(imt,jmt,km)::akmmon,aktmon,aksmon
!lhl1204
!
      real(r4),dimension(jmt_global,km+1,2):: psi
      real(r4),dimension(imt,jmt_global):: bsf
      real(r4),dimension(jmt_global,2,NTRA):: mth,mth_adv,mth_dif,mth_adv_iso
!
      real(r4),dimension(imt,jmt,km,NTRA)::trendmon
      real(r4),dimension(imt,jmt,km,NTRA)::axmon,aymon,azmon
      real(r4),dimension(imt,jmt,km,NTRA)::dxmon,dymon,dzmon
      real(r4),dimension(imt,jmt,km)::penmon
!
      real(r4),dimension(imt,jmt,km,NTRA)::ddymon

#ifdef ISO
      real(r4),dimension(imt,jmt,km,NTRA)::axmon_iso,aymon_iso,azmon_iso
      real(r4),dimension(imt,jmt,km,NTRA)::dxmon_iso,dymon_iso,dzmon_iso
      real(r4),dimension(imt,jmt,km,NTRA)::aaymon_iso,ddymon_iso
#endif
!
!#ifdef SPMD
!      real(r4),dimension(imt,jmt_global)::z0mon_io,himon_io,hdmon_io
!!
!!linpf091126      
!      real(r4),dimension(imt,jmt_global)::lthfmon_io,sshfmon_io,lwvmon_io,swvmon_io
!      real(r4),dimension(imt,jmt_global)::sumon_io,svmon_io
!!linpf091126
!!
!      real(r4),dimension(imt,jmt_global,km)::wsmon_io,tsmon_io,ssmon_io,usmon_io,vsmon_io
!      real(r4),dimension(imt,jmt_global,2):: icmon_io
!#if (defined SMAG_OUT)
!      real(r4),dimension(imt,jmt_global,km)::a3mon_io
!#endif
!      real(r4),dimension(imt,jmt_global,NTRA)::netmon_io
!      real(r4),dimension(imt,jmt_global,km,NTRA)::trendmon_io
!      real(r4),dimension(imt,jmt_global,km,NTRA)::axmon_io,aymon_io,azmon_io
!      real(r4),dimension(imt,jmt_global,km,NTRA)::dxmon_io,dymon_io,dzmon_io
!      real(r4),dimension(imt,jmt_global,km)::penmon_io
!      real(r4),dimension(imt,jmt_global,km,NTRA)::ddymon_io
!#ifdef ISO
!      real(r4),dimension(imt,jmt_global,km,NTRA)::axmon_iso_io,aymon_iso_io,azmon_iso_io
!      real(r4),dimension(imt,jmt_global,km,NTRA)::dxmon_iso_io,dymon_iso_io,dzmon_iso_io
!      real(r4),dimension(imt,jmt_global,km,NTRA)::aaymon_iso_io,ddymon_iso_io
!#endif
!      real(r4),dimension(imt,jmt_global)::mldmon_io
!      real(r4),dimension(imt,jmt_global,km)::akmmon_io,aktmon_io,aksmon_io
!#endif
#if (defined SMAG_OUT)
      real(r4),dimension(imt,jmt,km)::a3mon
#endif
      real(r4), parameter :: spval =1.0e+35
end module output_mod
