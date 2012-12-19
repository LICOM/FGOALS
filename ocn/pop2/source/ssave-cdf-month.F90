!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine SSAVEMON
!     =================
!     output in NETcdf format
!     written by liu hai long 2001 jun
!     =================
!     output history (netcdf) and restart (binary) files
!     remove pre-compilation NETCDF/ALL/NORMAL, remove yearly mean part 
!     keep the diagnostic part, SPMD
!     written by liu hailong & Lin Pengfei 2012 July
!     =================
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use output_mod
use dyn_mod
use tracer_mod
use cdf_mod
use diag_mod
#ifdef SPMD
use msg_mod
#endif

      IMPLICIT none
#include <netcdf.inc>
 
      logical :: hist_output,rest_output
      CHARACTER ( LEN =   4 ) :: ftail
      CHARACTER ( LEN =  24 ) :: fname
      CHARACTER ( LEN =  15 ) :: fname1
      CHARACTER ( LEN =   8 ) :: dd
      CHARACTER ( LEN =   10 ) :: tt
      CHARACTER ( LEN =   5 ) :: zz
      INTEGER(r4)             :: vv(8)
      INTEGER :: nwmf
 
!---------------------------------------------------------------------
!     output monthly results
!---------------------------------------------------------------------
!    file name
 
      nwmf = iyfm
      write (ftail,'(i4.4)') nwmf
      fname1(1:5)='MMEAN'
      fname1(6:9)=ftail
      fname1(10:10)='-'
      write(fname1(11:12),'(i2.2)')mon0
      fname1(13:15)='.nc'
 
!      if (iday==imd) then
      if (mod (iyfm,io_hist)==0 ) then
         hist_output=.true.
      else
         hist_output=.false.
      endif
!
      if (mod ((month-1),io_rest)==0 ) then
         rest_output=.true.
      else
         rest_output=.false.
      endif

       !write(*,*) 'hist_output,io_hist',hist_output,io_hist
       if(mytid==0) write(*,*) 'ok inssavemon'

#ifdef SPMD
      if (hist_output) then
!      call msf(vsmon_io)
!      call barosf(usmon_io)
!      call diag_tracer(1)
!      call diag_tracer(2)
!      call diag_heat_transport(1)
!      call diag_heat_transport(2)
#else
      IF (hist_output) THEN
!      call msf(vsmon)
!      call barosf(usmon)
!      call diag_tracer(1)
!      call diag_tracer(2)
!      call diag_heat_transport(1)
!      call diag_heat_transport(2)
#endif 

#if (defined NETCDF) || (defined ALL)
!--------------------------------------------------------------
!     cdf output
!--------------------------------------------------------------
!     file defination
! enter define mode
       if (mytid==0) then
         iret = nf_create (fname1, NF_CLOBBER, ncid)
         CALL check_err (iret)
          write(*,*)'ok generate,ncid=',ncid
! define dimensions
         iret = nf_def_dim (ncid, 'lat', lat_len, lat_dim)
         CALL check_err (iret)
         iret = nf_def_dim (ncid, 'lon', lon_len, lon_dim)
         CALL check_err (iret)
 
         IF (mon0 == 12) THEN
            iret = nf_def_dim (ncid, 'lev', lev_len, lev_dim)
            CALL check_err (iret)
         ELSE
            iret = nf_def_dim (ncid, 'lev', klv, lev_dim)
            CALL check_err (iret)
         END IF
 
         iret = nf_def_dim (ncid, 'lev1', lev1_len, lev1_dim)
         CALL check_err (iret)
 
         iret = nf_def_dim (ncid, 'time', NF_UNLIMITED, time_dim)
!      iret = nf_def_dim(ncid, 'time', time_len, time_dim)
         CALL check_err (iret)
! define variables
         lat_dims (1) = lat_dim
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, lat_dims, lat_id)
         CALL check_err (iret)
!
         lon_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, lon_dims, lon_id)
         CALL check_err (iret)
!
         lev_dims (1) = lev_dim
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, lev_dims, lev_id)
         CALL check_err (iret)
!
         lev1_dims (1) = lev1_dim
         iret = nf_def_var (ncid, 'lev1', NF_REAL, lev1_rank, lev1_dims, lev1_id)
         CALL check_err (iret)
!
         time_dims (1) = time_dim
         iret = nf_def_var (ncid, 'time', NF_DOUBLE, time_rank, time_dims, time_id)
         CALL check_err (iret)
 
         z0_dims (3) = time_dim
         z0_dims (2) = lat_dim
         z0_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'z0', NF_REAL, z0_rank, z0_dims, z0_id)
         CALL check_err (iret)
 
         hi_dims (3) = time_dim
         hi_dims (2) = lat_dim
         hi_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'hi', NF_REAL, hi_rank, hi_dims, hi_id)
         CALL check_err (iret)
 
         hd_dims (3) = time_dim
         hd_dims (2) = lat_dim
         hd_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'hd', NF_REAL, hd_rank, hd_dims, hd_id)
         CALL check_err (iret)
 
         ic1_dims (3) = time_dim
         ic1_dims (2) = lat_dim
         ic1_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ic1', NF_REAL, ic1_rank, ic1_dims, ic1_id)
         CALL check_err (iret)
 
         ic2_dims (3) = time_dim
         ic2_dims (2) = lat_dim
         ic2_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ic2', NF_REAL, ic2_rank, ic2_dims, ic2_id)
         CALL check_err (iret)

         net1_dims (3) = time_dim
         net1_dims (2) = lat_dim
         net1_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'net1', NF_REAL, net1_rank, net1_dims, net1_id)
         CALL check_err (iret)
 
         net2_dims (3) = time_dim
         net2_dims (2) = lat_dim
         net2_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'net2', NF_REAL, net2_rank, net2_dims, net2_id)
         CALL check_err (iret)
!
         mld_dims (3) = time_dim
         mld_dims (2) = lat_dim
         mld_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'mld', NF_REAL, mld_rank, mld_dims, mld_id)
         CALL check_err (iret)
         akm_dims (4) = time_dim
         akm_dims (3) = lev_dim
         akm_dims (2) = lat_dim
         akm_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'akm', NF_REAL, akm_rank, akm_dims, akm_id)
         CALL check_err (iret)
         akt_dims (4) = time_dim
         akt_dims (3) = lev_dim
         akt_dims (2) = lat_dim
         akt_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'akt', NF_REAL, akt_rank, akt_dims, akt_id)
         CALL check_err (iret)
         aks_dims (4) = time_dim
         aks_dims (3) = lev_dim
         aks_dims (2) = lat_dim
         aks_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'aks', NF_REAL, aks_rank, aks_dims, aks_id)
         CALL check_err (iret)
         ts_dims (4) = time_dim
         ts_dims (3) = lev_dim
         ts_dims (2) = lat_dim
         ts_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ts', NF_REAL, ts_rank, ts_dims, ts_id)
         CALL check_err (iret)
         ss_dims (4) = time_dim
         ss_dims (3) = lev_dim
         ss_dims (2) = lat_dim
         ss_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ss', NF_REAL, ss_rank, ss_dims, ss_id)
         CALL check_err (iret)
         us_dims (4) = time_dim
         us_dims (3) = lev_dim
         us_dims (2) = lat_dim
         us_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'us', NF_REAL, us_rank, us_dims, us_id)
         CALL check_err (iret)
         vs_dims (4) = time_dim
         vs_dims (3) = lev_dim
         vs_dims (2) = lat_dim
         vs_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'vs', NF_REAL, vs_rank, vs_dims, vs_id)
         CALL check_err (iret)
         ws_dims (4) = time_dim
         ws_dims (3) = lev_dim
         ws_dims (2) = lat_dim
         ws_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ws', NF_REAL, ws_rank, ws_dims, ws_id)
         CALL check_err (iret)
         su_dims (3) = time_dim
         su_dims (2) = lat_dim
         su_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'su', NF_REAL, su_rank, su_dims, su_id)
         CALL check_err (iret) !Uwindstress

         sv_dims (3) = time_dim
         sv_dims (2) = lat_dim
         sv_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'sv', NF_REAL, sv_rank, sv_dims, sv_id)
         CALL check_err (iret) !Vwindstress
         lthf_dims (3) = time_dim
         lthf_dims (2) = lat_dim
         lthf_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'lthf', NF_REAL, lthf_rank, lthf_dims, lthf_id)
         CALL check_err (iret) !latent heat flux 
         sshf_dims (3) = time_dim
         sshf_dims (2) = lat_dim
         sshf_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'sshf', NF_REAL, sshf_rank, sshf_dims, sshf_id)
         CALL check_err (iret) !sensible heat flux 
         lwv_dims (3) = time_dim
         lwv_dims (2) = lat_dim
         lwv_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'lwv', NF_REAL, lwv_rank, lwv_dims, lwv_id)
         CALL check_err (iret) !longwave 
         swv_dims (3) = time_dim
         swv_dims (2) = lat_dim
         swv_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'swv', NF_REAL, swv_rank, swv_dims, swv_id)
         CALL check_err (iret) !shortwave 

         if (diag_msf) then
            psi_dims (3) = time_dim
            psi_dims (2) = lev1_dim
            psi_dims (1) = lat_dim
            iret = nf_def_var (ncid, 'psi', NF_REAL, psi_rank, psi_dims, psi_id)
            CALL check_err (iret)
         end if

         if (diag_bsf) then
            bsf_dims (3) = time_dim
            bsf_dims (2) = lat_dim
            bsf_dims (1) = lon_dim
            iret = nf_def_var (ncid, 'bsf', NF_REAL, bsf_rank, bsf_dims, bsf_id)
            CALL check_err (iret)
         end if
 
         if (diag_mth) then
            mth_dims (3) = time_dim
            mth_dims (2) = lat_dim
            mth_dims (1) = lon_dim
            iret = nf_def_var (ncid, 'mth', NF_REAL, mth_rank, mth_dims, mth_id)
            CALL check_err (iret)
         end if

#if (defined SMAG_OUT)
         am_dims (4) = time_dim
         am_dims (3) = lev_dim
         am_dims (2) = lat_dim
         am_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'am', NF_REAL, am_rank, am_dims, am_id)
         CALL check_err (iret)
#endif
 
! assign attributes
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 18, 'depth (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev1_id, 'long_name', 18, 'depth (on V grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev1_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, time_id, 'long_name', 4, 'time')
         CALL check_err (iret)
!      iret = nf_put_att_text(ncid, time_id, 'units', 21, 'days since 1001-01-01')
         iret = nf_put_att_text (ncid, time_id, 'units', 23, 'months since 0001-01-01')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, z0_id, 'long_name', 18, 'sea surface height')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, z0_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, z0_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, hi_id, 'long_name', 16, 'thickness of ice')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, hi_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, hi_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, hd_id, 'long_name', 28, 'thickness of ice in one grid')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, hd_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, hd_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, ic1_id, 'long_name', 56, 'total of number of levels involved in convection per day')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ic1_id, 'units', 6, 'levels')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ic1_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, ic2_id, 'long_name', 35, 'number of levels ventilated per day')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ic2_id, 'units', 6, 'levels')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ic2_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, net1_id, 'long_name', 21, 'net surface heat flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, net1_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, net1_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, net2_id, 'long_name', 21, 'net surface salt flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, net2_id, 'units', 7, 'psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, net2_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, mld_id, 'long_name', 17, 'mixed layer depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, mld_id, 'units', 1, 'm')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, mld_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, akm_id, 'long_name', 28, 'turbulent vertical viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, akm_id, 'units', 5, 'm^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, akm_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, akt_id, 'long_name', 33, 'turbulent heat vertical viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, akt_id, 'units', 5, 'm^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, akt_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, aks_id, 'long_name', 33, 'turbulent salt vertical viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, aks_id, 'units', 5, 'm^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, aks_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, ts_id, 'long_name', 11, 'temperature')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ts_id, 'units', 10, 'centigrade')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ts_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, ss_id, 'long_name', 8, 'salinity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ss_id, 'units', 3, 'psu')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ss_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, us_id, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, us_id, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, us_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, vs_id, 'long_name', 18, 'meridional current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, vs_id, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, vs_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, ws_id, 'long_name', 16, 'vertical current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ws_id, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ws_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, su_id, 'long_name', 11, 'Uwindstress')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, su_id, 'units', 2, 'Pa')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, su_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, sv_id, 'long_name', 11, 'Vwindstress')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, sv_id, 'units', 2, 'Pa')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, sv_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, lthf_id, 'long_name', 16, 'latent heat flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lthf_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, lthf_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)


         iret = nf_put_att_text (ncid, sshf_id, 'long_name', 18, 'sensible heat flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, sshf_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, sshf_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, lwv_id, 'long_name', 8, 'Longwave')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lwv_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, lwv_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, swv_id, 'long_name', 9, 'Shortwave')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, swv_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, swv_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
!
         if (diag_msf) then
             iret = nf_put_att_text (ncid, psi_id, 'long_name', 27, 'Meridioanl Stream Function')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, psi_id, 'units', 8, 'Sverdrup')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, psi_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
         end if

         if (diag_bsf) then
             iret = nf_put_att_text (ncid, bsf_id, 'long_name', 26, 'Barotropic Stream Function')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, bsf_id, 'units', 8, 'Sverdrup')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, bsf_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
         end if
!
         if (diag_mth) then
             iret = nf_put_att_text (ncid, mth_id, 'long_name', 27, 'Meridional Tracer Transport')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, mth_id, 'units', 9, 'PW or m/s')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, mth_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
         end if
!
#if (defined SMAG_OUT)
         iret = nf_put_att_text (ncid, am_id, 'long_name', 20, 'horizontal viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, am_id, 'units', 6, 'm**2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, am_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
#endif
!   define global attribute
         CALL date_and_time (dd,tt,zz,vv)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'title', 4, 'test')
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'history', 20, tt //'  '//dd)
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
         CALL check_err (iret)
! leave define mode
         iret = nf_enddef (ncid)
         CALL check_err (iret)

!----------------------------------------------------------
!     prepare data for storing
!----------------------------------------------------------
 
         t0_cdf = month -1 !need to test 
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev1_id, lev1)
         CALL check_err (iret)
         start1 (1)= 1
         count1 (1)= time_len
         iret = nf_put_vara_double (ncid, time_id,start1,count1,t0_cdf)
         CALL check_err (iret)
 
         endif !mytid==0

! store variables
         start3 (1)= 1
         start3 (2)= 1
         start3 (3)= 1
         count3 (1)= lon_len
         count3 (2)= lat_len
         count3 (3)= time_len
 
         allocate(buffer_r4(imt_global,jmt_global))
#ifdef SPMD
        call local_to_global_4d(z0mon,buffer_r4,1,1,1)
        if(mytid==0) then         
         iret = nf_put_vara_real (ncid,z0_id,start3, count3, buffer_r4)
         CALL check_err (iret)
        endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= z0mon (i,j) / (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               ENDIF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO
         iret = nf_put_vara_real (ncid,z0_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif

#ifdef SPMD
         call local_to_global_4d(himon,buffer_r4,1,1,1)

         if(mytid==0) then         
         iret = nf_put_vara_real (ncid,hi_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= himon (i,j) / (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               ENDIF
            END DO
         END DO

         DO j = 1,jst_global
            DO i = 1,imt_global
            buffer_r4 (i,j,1)=spval            
            END DO
         END DO
         iret = nf_put_vara_real (ncid,hi_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif
 
#ifdef SPMD
         call local_to_global_4d(hdmon,buffer_r4,1,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,hd_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = 1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= hdmon (i,j)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO

         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,hd_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif

#ifdef SPMD
         call local_to_global_4d(icmon(:,:,1),buffer_r4,1,1,1)

         if(mytid==0) then
         iret = nf_put_vara_real (ncid,ic1_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0

#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+11,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= icmon (i,j,1)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO

         iret = nf_put_vara_real (ncid,ic1_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif

#ifdef SPMD
         call local_to_global_4d(icmon(:,:,2),buffer_r4,1,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,ic2_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= icmon (i,j,2)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO
         iret = nf_put_vara_real (ncid,ic2_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif
 
#ifdef SPMD
         call local_to_global_4d(netmon(:,:,1),buffer_r4,1,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,net1_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= netmon (i,j,1)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,net1_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif

#ifdef SPMD
         call local_to_global_4d(netmon(:,:,2),buffer_r4,1,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,net2_id,start3, count3, buffer_r4)
         CALL check_err (iret)
          endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= netmon (i,j,2)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO
         iret = nf_put_vara_real (ncid,net2_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif

#ifdef SPMD
         call local_to_global_4d(mldmon,buffer_r4,1,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,mld_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= mldmon (i,j)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO
         iret = nf_put_vara_real (ncid,mld_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif
!-----------------------------------------linpf 2012Jul27--------------
!3D output
         allocate(buffer3_r4(imt_global,jmt_global,km))
         
         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= klv
         count4 (4)= time_len

#ifdef SPMD
         call local_to_global_4d(akmmon,buffer3_r4,klv,1,0)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,akm_id,start4, count4, buffer3_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 1,klv
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (viv (i,j,k) > 0.5) THEN
                  buffer3_r4 (i,j,k,1)= akmmon (i,j,k)/ (nmonth (mon0))
               ELSE
                  buffer3_r4 (i,j,k,1)= spval
               END IF
            END DO
         END DO
         END DO
         DO k = 1,klv
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer3_r4 (i,j,1)= spval
            END DO
         END DO
         END DO
         iret = nf_put_vara_real (ncid,akm_id,start4, count4, buffer3_r4)
         CALL check_err (iret)
#endif

#ifdef SPMD
         call local_to_global_4d(aktmon,buffer3_r4,klv,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,akt_id,start4, count4, buffer3_r4)
         CALL check_err (iret)
         endif !mytid==0
#else    
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 1,klv
         DO j = 1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,k) > 0.5) THEN
                  buffer3_r4 (i,j,k,1)= aktmon (i,j,k)/ (nmonth (mon0))
               ELSE
                  buffer3_r4 (i,j,k,1)= spval
               END IF
            END DO
         END DO
         END DO
         DO k = 1,klv
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer3_r4 (i,j,k,1)= spval
            END DO
         END DO
         END DO
         iret = nf_put_vara_real (ncid,akt_id,start4, count4, buffer3_r4)
         CALL check_err (iret)
#endif 
!AKS
#ifdef SPMD
         call local_to_global_4d(aksmon,buffer3_r4,klv,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,aks_id,start4, count4, buffer3_r4)
         CALL check_err (iret)
         endif !mytid==0
#else    
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 1,klv
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,k) > 0.5) THEN
                  buffer3_r4 (i,j,k,1)= aksmon (i,j,k)/ (nmonth (mon0))
               ELSE
                  buffer3_r4 (i,j,k,1)= spval
               END IF
            END DO
         END DO
         END DO
         DO k = 1,klv
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer3_r4 (i,j,k)= spval
            END DO
         END DO
         END DO
         iret = nf_put_vara_real (ncid,aks_id,start4, count4, buffer3_r4)
         CALL check_err (iret)
#endif 

#ifdef SPMD
         call local_to_global_4d(tsmon,buffer3_r4,klv,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,ts_id, start4, count4, buffer3_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 1,klv
            DO j = jst_global+1,jmt_global
               DO i = 1,imt_global
                  IF (vit (i,j,k) > 0.5) THEN
                     buffer3_r4 (i,j,k,1)= tsmon (i,j,k)/ (nmonth (mon0))
                  ELSE
                     buffer3_r4 (i,j,k,1)= spval
                  END IF
               END DO
            END DO
         END DO
         DO k = 1,klv
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer3_r4 (i,j,k,1)= spval
            END DO
         END DO
         END DO
 
         iret = nf_put_vara_real (ncid,ts_id, start4, count4, buffer3_r4)
         CALL check_err (iret)
#endif

!salinity 
#ifdef SPMD
         call local_to_global_4d(ssmon,buffer3_r4,klv,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,ss_id, start4, count4, buffer3_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 1,klv
            DO j = jst_global+1,jmt_global
               DO i = 1,imt_global
                  IF (vit (i,j,k) > 0.5) THEN
                     buffer3_r4 (i,j,k,1)= ssmon (i,j,k)*1000./ (nmonth (mon0)) +35.
                  ELSE
                     buffer3_r4 (i,j,k,1)= spval
                  END IF
               END DO
            END DO
         END DO
         DO k = 1,klv
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer3_r4 (i,j,k,1)= spval
            END DO
         END DO
         END DO
 
         iret = nf_put_vara_real (ncid,ss_id, start4, count4, buffer3_r4)
         CALL check_err (iret)
#endif
!vertical current
 
#ifdef SPMD
         call local_to_global_4d(wsmon,buffer3_r4,klv,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,ws_id, start4, count4, buffer3_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 1,klv
            DO j = jst_global+1,jmt_global
               DO i = 1,imt_global
                  IF (vit (i,j,k) > 0.5) THEN
                     buffer3_r4 (i,j,k,1)= wsmon (i,j,k)/ (nmonth (mon0))
                  ELSE
                     buffer3_r4 (i,j,k,1)= spval
                  END IF
               END DO
            END DO
         END DO
         DO k = 1,klv
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer3_r4 (i,j,k,1)= spval
            END DO
         END DO
         END DO
 
         iret = nf_put_vara_real (ncid,ws_id, start4, count4, buffer3_r4)
         CALL check_err (iret)
#endif
!zonal curret 
#ifdef SPMD
         call local_to_global_4d(usmon,buffer3_r4,klv,1,0)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,us_id, start4, count4, buffer3_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 1,klv
            DO j = jst_global+1,jmt_global
               DO i = 1,imt_global
                  IF (viv (i,j,k) > 0.5) THEN
                     buffer3_r4 (i,j,k,1)= usmon (i,j,k)/ (nmonth (mon0))
                  ELSE
                     buffer3_r4 (i,j,k,1)= spval
                  END IF
               END DO
            END DO
         END DO
         DO k = 1,klv
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer3_r4 (i,j,k,1)= spval
            END DO
         END DO
         END DO
 
         iret = nf_put_vara_real (ncid,us_id, start4, count4, buffer3_r4)
         CALL check_err (iret)
#endif

!meridional curret 
#ifdef SPMD
         call local_to_global_4d(vsmon,buffer3_r4,klv,1,0)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,vs_id, start4, count4, buffer3_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 1,klv
            DO j = jst_global+1,jmt_global
               DO i = 1,imt_global
                  IF (viv (i,j,k) > 0.5) THEN
                     buffer3_r4 (i,j,k,1)= - vsmon (i,j,k)/ (nmonth (mon0))
                  ELSE
                     buffer3_r4 (i,j,k,1)= spval
                  END IF
               END DO
            END DO
         END DO
         DO k = 1,klv
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer3_r4 (i,j,k,1)= spval
            END DO
         END DO
         END DO
 
         iret = nf_put_vara_real (ncid,vs_id, start4, count4, buffer3_r4)
         CALL check_err (iret)
#endif
!Taux 
#ifdef SPMD
         call local_to_global_4d(sumon,buffer_r4,1,1,0)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,su_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (viv (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= sumon (i,j)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO

         iret = nf_put_vara_real (ncid,su_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif
!Tauy 
#ifdef SPMD
         call local_to_global_4d(svmon,buffer_r4,1,1,0)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,sv_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (viv (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= svmon (i,j)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO

         iret = nf_put_vara_real (ncid,sv_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif
!for latent heat flux
#ifdef SPMD
         call local_to_global_4d(lthfmon,buffer_r4,1,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,lthf_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= lthfmon (i,j)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO

         iret = nf_put_vara_real (ncid,lthf_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif
!for sensible heat flux
#ifdef SPMD
         call local_to_global_4d(sshfmon,buffer_r4,1,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,sshf_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= sshfmon (i,j)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO

         iret = nf_put_vara_real (ncid,sshf_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif
!for long wave radiation
#ifdef SPMD
         call local_to_global_4d(lwvmon,buffer_r4,1,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,lwv_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= lwvmon (i,j)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO

         iret = nf_put_vara_real (ncid,lwv_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif
!for shortwave radiation
#ifdef SPMD
         call local_to_global_4d(swvmon,buffer_r4,1,1,1)
         if(mytid==0) then
         iret = nf_put_vara_real (ncid,swv_id,start3, count3, buffer_r4)
         CALL check_err (iret)
         endif !mytid==0
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = jst_global+1,jmt_global
            DO i = 1,imt_global
               IF (vit (i,j,1) > 0.5) THEN
                  buffer_r4 (i,j,1)= swvmon (i,j)/ (nmonth (mon0))
               ELSE
                  buffer_r4 (i,j,1)= spval
               END IF
            END DO
         END DO
         DO j = 1,jst_global
            DO i = 1,imt_global
                  buffer_r4 (i,j,1)= spval
            END DO
         END DO

         iret = nf_put_vara_real (ncid,swv_id,start3, count3, buffer_r4)
         CALL check_err (iret)
#endif
!linpf 20120728
         if(mytid==0) write(*,*)'before diag_msf'
!linpf 20120728
       if(mytid==0) then

       if (diag_msf) then
          do k=1,km+1
          do j=1,jmt_global
             if (psi(j,k,1)<10000.) then
                t2z_cdf(j,k,1)=psi(j,k,1)/(nmonth (mon0))
             else
                t2z_cdf(j,k,1)=spval
             end if
          end do
          end do
!!
          start3 (1)= 1
          start3 (2)= 1
          start3 (3)= 1
          count3 (1)= lat_len
          count3 (2)= km+1
          count3 (3)= time_len
          iret = nf_put_vara_real (ncid,psi_id, start3, count3, t2z_cdf)
          CALL check_err (iret)
       end if
!
       if (diag_bsf) then
!!
          do j=1,jmt_global
          do i=1,imt
#ifdef SPMD
!             if (viv_global(i,j,1)>0.5) then !linpf 2012Jul27
             if (viv(i,j,1)>0.5) then
#else
             if (viv(i,j,1)>0.5) then
#endif
                t2_cdf(i,j,1)=bsf(i,j)/(nmonth (mon0))
             else
                t2_cdf(i,j,1)=spval
             end if
          end do
          end do
!!
          start3 (1)= 1
          start3 (2)= 1
          start3 (3)= 1
          count3 (1)= lon_len
          count3 (2)= lat_len
          count3 (3)= 1
          iret = nf_put_vara_real (ncid,bsf_id, start3, count3, t2_cdf)
          CALL check_err (iret)
       end if !diag_bsf
!!
       if (diag_mth) then
!!         
          do i=1,2
          do j=1,jmt_global
                t1_cdf(i,j,1)=mth(j,i,1)
          end do
          end do
          do i=3,4
          do j=1,jmt_global
                t1_cdf(i,j,1)=mth(j,i-3+1,2)
          end do
          end do

          do i=5,6
          do j=1,jmt_global
                t1_cdf(i,j,1)=mth_adv(j,i-5+1,1)
          end do
          end do
          do i=7,8
          do j=1,jmt_global
                t1_cdf(i,j,1)=mth_adv(j,i-7+1,2)
          end do
          end do
!!
          do i=9,10
          do j=1,jmt_global
                t1_cdf(i,j,1)=mth_dif(j,i-9+1,1)
          end do
          end do
          do i=11,12
          do j=1,jmt_global
                t1_cdf(i,j,1)=mth_dif(j,i-11+1,2)
          end do
          end do
!
          do i=13,14
          do j=1,jmt_global
                t1_cdf(i,j,1)=mth_adv_iso(j,i-13+1,1)
          end do
          end do
           do i=15,16
          do j=1,jmt_global
                t1_cdf(i,j,1)=mth_adv_iso(j,i-15+1,2)
          end do
          end do
!
          start3 (1)= 1
          start3 (2)= 1
          start3 (3)= 1
          count3 (1)= lon_len
          count3 (2)= lat_len
          count3 (3)= 1
          iret = nf_put_vara_real (ncid,mth_id, start3, count3, t1_cdf)
          CALL check_err (iret)

        end if !diag_mth
       endif !mytid==0
 
#if (defined SMAG_OUT)
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = 1,klv
            DO j = 1,jmt_global
               DO i = 1,imt
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i,j,k,1)= am3mon_io (i,j,k)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i,j,k,1)= am3mon (i,j,k)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i,j,k,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,am_id, start4, count4, t3_cdf)
         CALL check_err (iret)
#endif
#endif

     if(mytid==0) then 
         iret = nf_CLOSE (ncid)
         CALL check_err (iret)
     endif

     END IF
 
 
         CALL mm00 (klv)
!lhl20120728      IF (mod ( (month -1),12) == 0)THEN
      IF (mod ( (month -1),io_rest) == 0)THEN
         CALL yy00
      END IF
!

!      if (rest_output) then
!#ifdef SPMD
!      if(mytid==0)then
!         fname(1:8)='fort.22.'
!         fname(9:12)=ftail
!         fname(13:13)='-'
!         if ( mon0 < 12) then
!             write(fname(14:15),'(i2.2)')mon0+1
!         else
!             write(fname(14:15),'(i2.2)')mon0-11
!             write(fname(9:12),'(i4.4)')nwmf+1
!         end if
!         fname(16:24)='-01-00000'
!         OPEN (90,file = fname,form ='unformatted',status ='unknown')
!         write (90) h0_io,u_io,v_io,at_io,hi_io,itice_io,alead_io,month
!      endif
!#else
!         fname(1:8)='fort.22.'
!         fname(9:12)=ftail
!         fname(13:13)='-'
!         if ( mon0 < 12) then
!             write(fname(14:15),'(i2.2)')mon0+1
!         else
!             write(fname(14:15),'(i2.2)')mon0-11
!             write(fname(9:12),'(i4.4)')nwmf+1
!         end if
!         fname(16:24)='-01-00000'
!         OPEN (90,file = fname,form ='unformatted',status ='unknown')
!         WRITE (90) h0,u,v,at,hi,itice,alead,month
!#endif
!!
!         if (mytid==0) then
!#ifdef COUP
!             write(90)t_cpl_io,s_cpl_io,u_cpl_io,v_cpl_io,dhdx_io,dhdy_io,q_io
!#endif
!             close (90)
!         end if
! 
!#ifdef SPMD
!      if(mytid==0)then
!         open(22,file='fort.22',form='unformatted')
!         rewind 22
!         write (22) h0_io,u_io,v_io,at_io,hi_io,itice_io,alead_io,month
!      endif
!#else
!      REWIND 22
!      WRITE (22) h0,u,v,at,hi,itice,alead,month
!#endif
!!
!      if (mytid==0) then
!#ifdef COUP
!         write(22)t_cpl_io,s_cpl_io,u_cpl_io,v_cpl_io,dhdx_io,dhdy_io,q_io
!#endif
!         close(22)
!      end if

!      end if
! 
!---------------------------------------------------------------------
!     reset some arrays
!---------------------------------------------------------------------
! 
!      CALL vinteg (u,ub)
!      CALL vinteg (v,vb)
!! 
!!$OMP PARALLEL DO PRIVATE (k,j,i)
!      DO k = 1,km
!         DO j = 1,jmt ! Dec. 5, 2002, Yongqiang Yu
!            DO i = 1,imt
!               up (i,j,k) = u (i,j,k)
!               vp (i,j,k) = v (i,j,k)
!               utf (i,j,k) = u (i,j,k)
!               vtf (i,j,k) = v (i,j,k)
!               atb (i,j,k,1) = at (i,j,k,1)
!               atb (i,j,k,2) = at (i,j,k,2)
!            END DO
!         END DO
!      END DO
! 
!!$OMP PARALLEL DO PRIVATE (j,i)
!      DO j = 1,jmt ! Dec. 5, 2002, Yongqiang Yu
!         DO i = 1,imt
!            h0p (i,j)= h0 (i,j)
!            ubp (i,j)= ub (i,j)
!            vbp (i,j)= vb (i,j)
!            h0f (i,j)= h0 (i,j)
!            h0bf (i,j)= h0 (i,j)
!         END DO
!      END DO
! 
      deallocate (buffer_r4,buffer3_r4)
      RETURN
      END 
