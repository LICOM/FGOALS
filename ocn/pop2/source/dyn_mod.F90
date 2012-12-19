!  CVS: $Id: dyn_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module dyn_mod
#include <def-undef.h>
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     U V T S H0 W RHO
!     ------------------------------------------------------------------
      real(r8),dimension(imt,jmt)::ub,vb,ubp,vbp,h0p
      real(r8),dimension(imt,jmt,km)::up,vp
      real(r8),dimension(imt,jmt,kmp1)::ws
      real(r8),dimension(imt,jmt)::h0l,h0f,h0bl,h0bf
      real(r8),dimension(imt,jmt,km)::utl,utf,vtl,vtf
!#ifdef SPMD
!      real(r8),dimension(imt_global,jmt_global,km)::u_io,v_io
!      real(r8),dimension(imt_global,jmt_global)::h0_io
!mohr
      real(r8),allocatable,dimension(:,:):: buffer
      real(r8),allocatable,dimension(:,:)::h0
      real(r8),allocatable,dimension(:,:,:)::u,v
!
#ifdef COUP
      real(r4),dimension(imt_global,jmt_global)::t_cpl_io,s_cpl_io,u_cpl_io,v_cpl_io,dhdx_io,dhdy_io,q_io
#endif
!#endif
!
!
!     ------------------------------------------------------------------
!     Pressure gradient
!     ------------------------------------------------------------------
      real(r8),dimension(:,:,:),allocatable::gg,dlu,dlv
      real(r8),dimension(:,:),allocatable::dlub,dlvb
end module dyn_mod

