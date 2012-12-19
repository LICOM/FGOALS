!  CVS: $Id: cdf_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
module cdf_mod
#include <def-undef.h>
use precision_mod
use param_mod
!     this head file includes variables related to Netcdf output.
!     written by liu hai long 2001, jun 
#if (defined NETCDF) || (defined ALL)
!     error status return
      integer:: iret

!     file id
      integer:: ncid

!     dimension id
      integer:: lon_dim,lat_dim,lev1_dim,lev_dim,time_dim

!     dimension lenth
      integer, parameter :: lon_len=imt_global,lat_len=jmt_global,lev_len=km,lev1_len=km+1,time_len=1

!     variable id
#if (defined SMAG_OUT)
      integer, parameter :: am_rank=4
      integer:: am_id,am_rank,am_dims(am_rank)
#endif
      integer::  lat_id,lon_id,lev1_id,lev_id,time_id,z0_id,hi_id,hd_id, &
               ic1_id,ic2_id, net1_id, net2_id, ts_id,  ss_id,us_id,vs_id,ws_id, psi_id, bsf_id,&
               mth_id,mld_id,akm_id,akt_id,aks_id &
               ,su_id,sv_id,sshf_id,lthf_id,lwv_id,swv_id&
             ,h0_id,u_id,v_id,t_id,s_id,h1_id,h2_id,h3_id,month_id, day_id

!     variable rank
!YU
      integer, parameter ::lat_rank=1,lon_rank=1,lev_rank=1,time_rank=1,z0_rank=3,hi_rank=3, &
                            hd_rank=3,ic1_rank=3,ic2_rank=3, net1_rank=3,net2_rank=3,ts_rank=4,ss_rank=4, &
                            us_rank=4, vs_rank=4, ws_rank=4, psi_rank=3,lev1_rank=1, bsf_rank=3,&
               mth_rank=3,mld_rank=3,akm_rank=4,akt_rank=4,aks_rank=4 &
              ,su_rank=3,sv_rank=3,sshf_rank=3,lthf_rank=3,lwv_rank=3,swv_rank=3&
             ,h0_rank=2,u_rank=3,v_rank=3,t_rank=3,s_rank=3,h1_rank=2,h2_rank=2,h3_rank=2,month_rank=1, day_rank=1
!
!     variable shapes
      integer::  lat_dims(lat_rank),lon_dims(lon_rank),lev_dims(lev_rank),time_dims(time_rank), &
                z0_dims( z0_rank), hi_dims( hi_rank), hd_dims( hd_rank), ic1_dims( ic1_rank), &
               ic2_dims(ic2_rank),net1_dims(net1_rank), net2_dims(net2_rank),ts_dims( ts_rank), ss_dims( ss_rank), &
                us_dims( us_rank), vs_dims( vs_rank), ws_dims( ws_rank),psi_dims(psi_rank), &
              lev1_dims(lev1_rank), bsf_dims(bsf_rank),&
              mth_dims(mth_rank),mld_dims(mld_rank),akm_dims(akm_rank),akt_dims(akt_rank),aks_dims(aks_rank), day_dims(day_rank), month_dims(month_rank) &
             ,su_dims(su_rank),sv_dims(sv_rank)&
             ,sshf_dims(sshf_rank),lthf_dims(lthf_rank),lwv_dims(lwv_rank),swv_dims(swv_rank)&
             ,h0_dims(h0_rank),u_dims(u_rank),v_dims(v_rank),t_dims(t_rank),s_dims(s_rank)&
             ,h1_dims(h1_rank),h2_dims(h2_rank),h3_dims(h3_rank)

      real(r8) t0_cdf 
      real(r4),dimension(imt_global,jmt_global,1):: t2_cdf 
      real(r4),dimension(imt_global,jmt_global,1)::  t2z_cdf,t1_cdf 
      real(r4),dimension(imt_global,jmt_global,km,1)::  t3_cdf !linpf 2012Jul27
      real(r4),allocatable,dimension(:,:):: buffer_r4 !linpf 2012Jul27
      real(r4),allocatable,dimension(:,:,:):: buffer3_r4 !linpf 2012Jul27
!     variables
!    start and count
      integer:: start1(1),count1(1)
      integer:: start2(2),count2(2)
      integer:: start3(3),count3(3)
      integer:: start4(4),count4(4)
#endif
end module cdf_mod
