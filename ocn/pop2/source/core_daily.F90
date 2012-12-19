!     =================
      SUBROUTINE CORE_DAILY(TNUM)
!     =================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use forc_mod
use dyn_mod, only: u,v,buffer
use tracer_mod, only: at
#if ( defined SPMD )
use msg_mod
#endif
      IMPLICIT NONE
#include <netcdf.inc>

      integer :: mon_day,irec,tnum
      real(r8),dimension(s_imt,s_jmt) :: t10,u10,v10,slp,q10,swhf,lwhf
      real(r8),dimension(s_imt,s_jmt) :: precr,precs

      real(r8),dimension(imt,jmt) :: model_sst,es,qs,zz,uu,vv,windx,windy,theta
      real(r8),dimension(imt,jmt) :: core_sensible,core_latent,core_tau
!
      real(r8), parameter :: tok=273.15
      real(r8), parameter :: epsln=1e-25
      real(r8),dimension(imt,jmt) :: tmp1,tmp2
      real(r8),dimension(imt_global,jmt_global) :: tmp3


      if ( TNUM.gt.1 ) goto 111

! decide recode number
!      IYFM,MON0,IDAY
      mon_day=0
      do i=1,mon0-1
      mon_day=mon_day+nmonth(i)
      enddo
!
! start from a 1000-year spinup
      irec=(iyfm-1)*365+mon_day+iday-365*(13-1)
!climatology forcing
!      irec=mon_day+iday

      if (mytid.eq.0) then
!      write(*,*) "iyfm=",iyfm
!      write(*,*) "mon0=",mon0
!      write(*,*) "iday=",iday
!      write(*,*) "mon_day=",mon_day
      write(*,*) "irec=",irec
!      write(*,*) s_imt,s_jmt
      endif

! read in core data
! note the dimensions are s_imt, s_jmt
#ifdef SPMD
      if (mytid.eq.0) then

      call read_core(irec,"t_10.db.1948-2007.daymean.05APR2010.nc",t10)
      call read_core(irec,"u_10.db.1948-2007.daymean.05APR2010.nc",u10)
      call read_core(irec,"v_10.db.1948-2007.daymean.05APR2010.nc",v10)
      call read_core(irec,"slp.db.1948-2007.daymean.05APR2010.nc",slp)
      call read_core(irec,"q_10.db.1948-2007.daymean.05APR2010.nc",q10)
      call read_core(irec,"swdn.db.1948-2007.daymean.05APR2010.nc",swhf)
      call read_core(irec,"lwdn.db.1948-2007.daymean.05APR2010.nc",lwhf)
      call read_core(irec,"rain.db.1948-2007.daymean.05APR2010.nc",precr)
      call read_core(irec,"snow.db.1948-2007.daymean.05APR2010.nc",precs)

!      call read_core(irec,"t_10.db.clim.daymean.15JUNE2009.nc",t10)
!      call read_core(irec,"u_10.db.clim.daymean.15JUNE2009.nc",u10)
!      call read_core(irec,"v_10.db.clim.daymean.15JUNE2009.nc",v10)
!      call read_core(irec,"slp.db.clim.daymean.15JUNE2009.nc",slp)
!      call read_core(irec,"q_10.db.clim.daymean.15JUNE2009.nc",q10)
!      call read_core(irec,"swdn.db.clim.daymean.15JUNE2009.nc",swhf)
!      call read_core(irec,"lwdn.db.clim.daymean.15JUNE2009.nc",lwhf)
!      call read_core(irec,"rain.db.clim.daymean.15JUNE2009.nc",precr)
!      call read_core(irec,"snow.db.clim.daymean.15JUNE2009.nc",precs)


      endif

! interplate to T grid
! the dimensions are imt and jmt_global
      if (mytid.eq.0) then
      call interplation(t10,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,tsa3(1,1,1))

      if (mytid.eq.0) then
      call interplation(u10,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,wspdu3(1,1,1))

      if (mytid.eq.0) then
      call interplation(v10,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,wspdv3(1,1,1))

      if (mytid.eq.0) then
      call interplation(slp,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,psa3(1,1,1))

      if (mytid.eq.0) then
      call interplation(q10,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,qar3(1,1,1))

      if (mytid.eq.0) then
      call interplation(swhf,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,swv3(1,1,1))

      if (mytid.eq.0) then
      call interplation(lwhf,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,lwv3(1,1,1))

      if (mytid.eq.0) then
      call interplation(precr,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,rain3(1,1,1))

      if (mytid.eq.0) then
      call interplation(precs,tmp3)
!      do j=1,jmt_global
!      do i=1,imt_global
!         tmp3(i,j)=tmp3(i,j)*vit_global_surface(i,j)
!      end do
!      end do
      endif
      call global_distribute(tmp3,snow3(1,1,1))


!      call interplation(t10,tsa3_io(1,1,1))
!      call interplation(u10,wspdu3_io(1,1,1))
!      call interplation(v10,wspdv3_io(1,1,1))
!      call interplation(slp,psa3_io(1,1,1))
!      call interplation(q10,qar3_io(1,1,1))
!      call interplation(swhf,swv3_io(1,1,1))
!      call interplation(lwhf,lwv3_io(1,1,1))
!      call interplation(precr,rain3_io(1,1,1))
!      call interplation(precs,snow3_io(1,1,1))
!
!      do j=1,jmt_global
!         do i=1,imt
!            tsa3_io(i,j,1)=tsa3_io(i,j,1)*vit_global(i,j,1)
!            wspdu3_io(i,j,1)=wspdu3_io(i,j,1)*vit_global(i,j,1)
!            wspdv3_io(i,j,1)=wspdv3_io(i,j,1)*vit_global(i,j,1)
!            psa3_io(i,j,1)=psa3_io(i,j,1)*vit_global(i,j,1)
!            qar3_io(i,j,1)=qar3_io(i,j,1)*vit_global(i,j,1)
!            swv3_io(i,j,1)=swv3_io(i,j,1)*vit_global(i,j,1)
!            lwv3_io(i,j,1)=lwv3_io(i,j,1)*vit_global(i,j,1)
!            rain3_io(i,j,1)=rain3_io(i,j,1)*vit_global(i,j,1)
!            snow3_io(i,j,1)=snow3_io(i,j,1)*vit_global(i,j,1)
!         end do
!      end do
!
!      endif
!
!       call global_to_local_4d(tsa3_io(1,1,1),tsa3(1,1,1),1,1)
!       call global_to_local_4d(wspdu3_io(1,1,1),wspdu3(1,1,1),1,1)
!       call global_to_local_4d(wspdv3_io(1,1,1),wspdv3(1,1,1),1,1)
!       call global_to_local_4d(psa3_io(1,1,1),psa3(1,1,1),1,1)
!       call global_to_local_4d(qar3_io(1,1,1),qar3(1,1,1),1,1)
!       call global_to_local_4d(swv3_io(1,1,1),swv3(1,1,1),1,1)
!       call global_to_local_4d(lwv3_io(1,1,1),lwv3(1,1,1),1,1)
!       call global_to_local_4d(rain3_io(1,1,1),rain3(1,1,1),1,1)
!       call global_to_local_4d(snow3_io(1,1,1),snow3(1,1,1),1,1)

#endif

111    continue

         do j = jsm,jem
            do i = 2,imm
              uu(i,j)=vit(i,j,1)*(u(i,j,1)+u(i+1,j,1)+u(i,j-1,1)+u(i+1,j-1,1)) &
              /(viv(i,j,1)+viv(i+1,j,1)+viv(i,j-1,1)+viv(i+1,j-1,1)+epsln)
              vv(i,j)=vit(i,j,1)*(v(i,j,1)+v(i+1,j,1)+v(i,j-1,1)+v(i+1,j-1,1)) &
              /(viv(i,j,1)+viv(i+1,j,1)+viv(i,j-1,1)+viv(i+1,j-1,1)+epsln)
            end do
         end do

     if (nx_proc == 1) then
            do j=jst,jem
              uu(1,j)=uu(imm,j)
              uu(imt,j)=uu(2,j)
              vv(1,j)=vv(imm,j)
              vv(imt,j)=vv(2,j)
            end do
     endif

#ifdef SPMD
       call exch_boundary(uu,1)
       call exch_boundary(vv,1)
#endif

! transfer core data to what the subroutine need
      do j=1,jmt
         do i=1,imt
! relative speed to surface currents
         windx(i,j)=(wspdu3(i,j,1)-uu(i,j))*vit(i,j,1)
         windy(i,j)=(wspdv3(i,j,1)-vv(i,j))*vit(i,j,1)
!  1.0 is from mom4
         wspd3(i,j,1)=sqrt(windx(i,j)**2+windy(i,j)**2+1.0)*vit(i,j,1)
! using a transient temperature, not daily mean
         model_sst(i,j)=(at(i,j,1,1)+tok)*vit(i,j,1)
         zz(i,j)=10.
         qs(i,j)=0.98*640380*exp(-5107.4/model_sst(i,j))/1.22*vit(i,j,1)
! temperature to potential temperature
         theta(i,j)=tsa3(i,j,1)*(100000.0/psa3(i,j,1))**0.286*vit(i,j,1)
         end do
      end do


! compute heat flux
       call ncar_ocean_fluxes(wspd3(1,1,1),theta(1,1),model_sst(1,1),qar3(1,1,1),qs(1,1),zz(1,1),vit(1,1,1),&
           core_sensible(1,1),core_latent(1,1),core_tau,ustar(1,1))

       do j=1,jmt
          do i=1,imt
            sshf(i,j)=core_sensible(i,j)*vit(i,j,1)*(1.0d0-seaice(i,j))
            lthf(i,j)=(core_latent(i,j)-snow3(i,j,1)*3.335e+5)*vit(i,j,1)*(1.0d0-seaice(i,j))
            lwv(i,j)=(0.95*lwv3(i,j,1)-0.95*5.67E-8*model_sst(i,j)**4)*vit(i,j,1)*(1.0d0-seaice(i,j))
            tmp1(i,j)=core_tau(i,j)*windx(i,j)*vit(i,j,1)
            tmp2(i,j)=-core_tau(i,j)*windy(i,j)*vit(i,j,1)
            nswv(i,j)=(lwv(i,j)+sshf(i,j)+lthf(i,j))
            swv(i,j)=(1-0.066)*swv3(i,j,1)*vit(i,j,1)*(1.0d0-seaice(i,j))
            fresh(i,j)=-(core_latent(i,j)/(2.5e+6)+rain3(i,j,1)+snow3(i,j,1)+runoff(i,j))*vit(i,j,1)*(1.0d0-seaice(i,j))
            ustar(i,j)=ustar(i,j)*vit(i,j,1)
         end do
      end do

! tau to U/V grid
        do j = jst,jem
           do i = 2,imm
             su(i,j)= 0.25*(tmp1(i,j)+tmp1(i-1,j)+tmp1(i,j+1)+tmp1(i-1,j+1))*viv(i,j,1)
             sv(i,j)= 0.25*(tmp2(i,j)+tmp2(i-1,j)+tmp2(i,j+1)+tmp2(i-1,j+1))*viv(i,j,1)
           end do
        end do

     if (nx_proc == 1) then
            do j=jst,jem
           su(1,j)=su(imm,j)
           su(imt,j)=su(2,j)
           sv(1,j)=sv(imm,j)
           sv(imt,j)=sv(2,j)
            end do
     end if
#ifdef SPMD
       call exch_boundary(su,1)
       call exch_boundary(sv,1)
#endif
!       allocate(buffer(imt_global,jmt_global))
!       call local_to_global_4d_double(su,buffer,1,1)
!       if (mytid.eq.0) then
!       write(133,*) buffer
!       close(133)
!       endif
!       call local_to_global_4d_double(sv,buffer,1,1)
!       if (mytid.eq.0) then
!       write(134,*) buffer
!       close(134)
!       endif
!       deallocate(buffer)
!       stop


      return
      end subroutine CORE_DAILY


!---------------------------------------------
      subroutine read_core(nnn,fname,var)
!---------------------------------------------
use precision_mod
      use param_mod, only: s_imt,s_jmt
      implicit none
#include <netcdf.inc>

      integer :: start(3),count(3)
      integer :: ncid,iret,nnn
      real(r8) :: var(s_imt,s_jmt)
      character (len=180) :: fname

      start(1)=1;count(1)=s_imt
      start(2)=1;count(2)=s_jmt
      start(3)=nnn;count(3)=1

      iret=nf_open(fname,nf_nowrite,ncid)
      call check_err (iret)

      iret=nf_get_vara_double(ncid,1,start,count,var)
      call check_err (iret)

      iret=nf_close(ncid)
      call check_err (iret)

      return
      end subroutine read_core

!---------------------------------------------
      subroutine read_core1(nnn,fname,var)
!---------------------------------------------
use precision_mod
      use param_mod, only: s_imt,s_jmt
      implicit none
#include <netcdf.inc>

      integer :: start(3),count(3)
      integer :: ncid,iret,nnn
      real(r8) :: var(s_imt,s_jmt)
      character (len=180) :: fname

      start(1)=1;count(1)=s_imt
      start(2)=1;count(2)=s_jmt
      start(3)=nnn;count(3)=1

      iret=nf_open(fname,nf_nowrite,ncid)
      call check_err (iret)

      iret=nf_get_vara_double(ncid,3,start,count,var)
      call check_err (iret)

      iret=nf_close(ncid)
      call check_err (iret)


      return
      end subroutine read_core1

!---------------------------------------------
      subroutine interplation(source,object)
!---------------------------------------------
use precision_mod
      use param_mod, only: imt,jmt,imt_global,jmt_global,s_imt,s_jmt,mytid
      use pconst_mod, only: s_lon,s_lat
      use output_mod, only: spval
      implicit none

      integer, parameter :: iwk=s_imt+2,jwk=s_jmt+2
      real(r8) :: source(s_imt,s_jmt)
      real(r8) :: object(imt_global,jmt_global)
      real(r8) :: s_work(iwk,jwk),s_wx(iwk),s_wy(jwk)
      integer :: i,j

       do j=1,s_jmt
       do i=1,s_imt
       if(source(i,j).lt.-1.0e+30)then
        source(i,j)=spval
       endif
       enddo
       enddo

       do i=2,iwk-1
       s_wx(i)=s_lon(i-1)
       enddo
       s_wx(1)=s_wx(2)-1.875
       s_wx(iwk)=s_wx(iwk-1)+1.875

! s_lat is already from north to south!
       do j=2,jwk-1
       s_wy(j)=s_lat(j-1)
       enddo
       s_wy(1)=s_wy(2)+1.875
       s_wy(jwk)=s_wy(jwk-1)-1.875
!	if(mytid.eq.0) write(*,*)s_wx
!	if(mytid.eq.0) write(*,*)s_wy
!	stop
!
! source from north to south!
      do j=1,s_jmt
      do i=1,s_imt
      s_work(i+1,j+1)=source(i,s_jmt-j+1)
      enddo
!wrong      s_work(1,j+1)=source(imt,s_jmt-j+1)
      s_work(1,j+1)=source(s_imt,s_jmt-j+1)
      s_work(iwk,j+1)=source(1,s_jmt-j+1)
      enddo
!
      do i=1,s_imt
      s_work(i+1,1)=source(i,2)
      s_work(i+1,jwk)=source(i,1)
      enddo
      s_work(  1,  1)=source(s_imt,  2)
      s_work(iwk,jwk)=source(  1,1)
!
      call hsetup(s_work,iwk,jwk,s_wx,s_wy,object(1,1))

      return
      end subroutine interplation

!---------------------------------------------
      subroutine hsetup(a,mx,my,alon,alat,b)
!---------------------------------------------
!     input : a
!     output: b
!     A bi-linear interpolation will be used for the initial guess, then
!     a refilling procedure be called to redefine the missing data.
!     mx,my            x and y grids number of a
!     alon,alat        lontitude and latitude of a
!     spval            missing flag
!     nf               1 for T grid; 0 for U grid
!
use precision_mod
      use param_mod, only: imt_global,jmt_global,s_imt,s_jmt
      use pconst_mod, only: lon,lat,vit_global_surface
      use output_mod, only: spval
      implicit none
!
      integer :: mx,my,ic,jc,ip,jp,i,j
!      integer :: ind(imt_global,jmt_global)
      real(r8):: b(imt_global,jmt_global)
      real(r8):: a(mx,my),alat(my),alon(mx)
      real(r8), parameter :: isp=99999.0d0
      real(r8):: r1,r2,s1,s2,b1,b2
!
!   Initiallize
!
!  ---nf = 1 on T girds

!  ---b results
      do j=1,jmt_global
      do i=1,imt_global
      b(i,j)=spval
      enddo
      enddo
!
      do 100 j=1,jmt_global
      do 100 i=1,imt_global-2
!
      if(vit_global_surface(i,j).lt.0.01) goto 100
!
!  ---find adjacent two grids on x direction (ic)
      ic=isp
      jc=isp
      do 45 ip=2,mx
      if(alon(ip-1).le.lon(i).and.alon(ip).ge.lon(i))then
      ic=ip
      goto 33
      endif
   45 continue
   33 continue

!
!  ---find adjacent two grids on y direction (jc)
      do 50 jp=2,my
      if(alat(jp).le.lat(j).and.alat(jp-1).ge.lat(j))then
      jc=jp
      goto 44
      endif
   50 continue
   44 continue

!
!  ---break if adjacent grids has no found
      if(ic.eq.isp.or.jc.eq.isp) then
      write(*,*)ic,jc
      write(*,*)i,j,alon(i),alat(j)
      write(*,*)lat
      stop 999
      endif
!
!  Bilinear interpolater
!
      r1=(lon(i)-alon(ic-1))/(alon(ic)-alon(ic-1))
      r2=(alon(ic)-lon(i))/(alon(ic)-alon(ic-1))
      s1=(lat(j)-alat(jc-1))/(alat(jc)-alat(jc-1))
      s2=(alat(jc)-lat(j))/(alat(jc)-alat(jc-1))
!


      if(a(ic-1,jc-1).eq.spval.or.a(ic,jc-1).eq.spval.or.&
       a(ic-1,jc).eq.spval.or.a(ic,jc).eq.spval)then

      b1=spval
      b2=spval
      b(i,j)=spval
      else

      b1=a(ic-1,jc-1)*r2+a(ic,jc-1)*r1
      b2=a(ic-1,jc)*r2+a(ic,jc)*r1

      b(i,j)=s1*b2+s2*b1

      endif
!
  100 continue
!
!    set cyclic b. c.
!
      do j=1,jmt_global
      b(imt_global-1,j)     = b(1,j)
      b(imt_global,j)       = b(2,j)
      enddo
!
      return
      end subroutine hsetup

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Over-ocean fluxes following Large and Yeager (used in NCAR models)           !
! Coded by Mike Winton (Michael.Winton@noaa.gov) in 2004
!
! A bug was found by Laurent Brodeau (brodeau@gmail.com) in 2007.
! Stephen.Griffies@noaa.gov updated the code with the bug fix. 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
subroutine ncar_ocean_fluxes (u_del, t, ts, q, qs, z, avail, &
                              sh,lh,tau,ustar)
!                              cd, ch, ce, ustar, bstar       )
use precision_mod
      use param_mod, only: imt,jmt,imt_global,jmt_global,s_imt,s_jmt,mytid
      implicit none
    real(r8)   , intent(in)   , dimension(imt,jmt) :: u_del, t, ts, q, qs, z
    real(r8)   , intent(in)   , dimension(imt,jmt) :: avail
    real(r8)   , intent(inout), dimension(imt,jmt) :: lh,sh,tau
    real(r8)   , dimension(imt,jmt) :: cd, ch, ce, ustar, bstar
!    real   , intent(inout), dimension(imt,jmt) :: cd, ch, ce, ustar, bstar

  real(r8) :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
  real(r8) :: cd_rt                                ! full drag coefficients @ z
  real(r8) :: zeta, x2, x, psi_m, psi_h            ! stability parameters
  real(r8) :: u, u10, tv, tstar, qstar, z0, xx, stab

  integer, parameter :: n_itts = 2
  real(r8), parameter :: grav = 9.80, vonkarm = 0.40,  L=2.5e6, cp=1000.5, r0=1.22
  integer               i, j, jj


  do j=1,jmt
  do i=1,imt
!  do i=1,size(u_del)
    if (avail(i,j) > 0.5 ) then
      tv = t(i,j)*(1+0.608*q(i,j));
      u = max(u_del(i,j), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
      u10 = u;                                                ! first guess 10m wind
    
      cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
      cd_n10_rt = sqrt(cd_n10);
      ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
      stab = 0.5 + sign(0.5,t(i,j)-ts(i,j))
      ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c
  
      cd(i,j) = cd_n10;                                         ! first guess for exchange coeff's at z
      ch(i,j) = ch_n10;
      ce(i,j) = ce_n10;
      do jj=1,n_itts                                           ! Monin-Obukhov iteration
        cd_rt = sqrt(cd(i,j));
        ustar(i,j) = cd_rt*u;                                   ! L-Y eqn. 7a
        tstar    = (ch(i,j)/cd_rt)*(t(i,j)-ts(i,j));                ! L-Y eqn. 7b
        qstar    = (ce(i,j)/cd_rt)*(q(i,j)-qs(i,j));                ! L-Y eqn. 7c
        bstar(i,j) = grav*(tstar/tv+qstar/(q(i,j)+1/0.608));

        zeta     = vonkarm*bstar(i,j)*z(i,j)/(ustar(i,j)*ustar(i,j)); ! L-Y eqn. 8a
        zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
        x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
        x2 = max(x2, 1.0);                                    ! undocumented NCAR
        x = sqrt(x2);
    
        if (zeta > 0) then
          psi_m = -5*zeta;                                    ! L-Y eqn. 8c
          psi_h = -5*zeta;                                    ! L-Y eqn. 8c
        else
          psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
          psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
        end if
    
        u10 = u/(1+cd_n10_rt*(log(z(i,j)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9


        cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
        cd_n10_rt = sqrt(cd_n10);
        ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
        stab = 0.5 + sign(0.5,zeta)
        ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
        z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic
    
        xx = (log(z(i,j)/10)-psi_m)/vonkarm;
        cd(i,j) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
        xx = (log(z(i,j)/10)-psi_h)/vonkarm;
!!$        ch(i,j) = ch_n10/(1+ch_n10*xx/cd_n10_rt)**2;                !       b (bug)
!!$        ce(i,j) = ce_n10/(1+ce_n10*xx/cd_n10_rt)**2;                !       c (bug)
        ch(i,j) = ch_n10/(1+ch_n10*xx/cd_n10_rt)*sqrt(cd(i,j)/cd_n10) ! 10b (corrected code aug2007)
        ce(i,j) = ce_n10/(1+ce_n10*xx/cd_n10_rt)*sqrt(cd(i,j)/cd_n10) ! 10c (corrected code aug2007)
      end do
    end if
  end do
  end do

    sh=r0*cp*ch*(t-ts)*u_del
    lh=r0*ce*l*(q-qs)*u_del
    tau=r0*cd*u_del

    return
end subroutine ncar_ocean_fluxes
