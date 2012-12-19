!  CVS: $Id: inirun.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE INIRUN
!     =================
!     INITIALIZING FOR ALL PHYSICAL FIELDS
 
#include <def-undef.h>
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
use shr_msg_mod
use shr_mpi_mod
use shr_sys_mod
use buf_mod
use control_mod


      IMPLICIT NONE

  !----- local  ------
  integer            :: fid    ! nc domain file ID
  integer            :: dimid  ! nc dimension id
  integer            :: vid    ! nc variable ID
  integer            :: rcode  ! nc return code
  integer            :: ntim   ! temporary
 
#include <netcdf.inc>
!
!     Define Variables.
      integer*4   :: ncid, iret
      integer*4,  dimension(4) :: start(4)
      integer*4,  dimension(4) :: count(4)
      character (len=18) :: fname

      INTEGER :: NMFF

      allocate(h0(imt,jmt),u(imt,jmt,km),v(imt,jmt,km),at(imt,jmt,km,ntra),hi(imt,jmt),itice(imt,jmt),alead(imt,jmt))
      allocate(buffer(imt_global,jmt_global),buffer_real4(imt_global,jmt_global))


      if (mytid==0)then
          write(6,*)"BEGINNING-----INIRUN !"
          open (17,file='rpointer.ocn',form='formatted')
          read(17,'(a18)') fname
          close(17)
      endif 
 
#ifdef COUP
!
     ncpl=1
!
  if (mytid==0 ) then
     call wrap_open ('domain_licom.nc', NF_NOWRITE, fid)

     write(6,*) 'read domain data...'
  call shr_sys_flush(6)

     ! obtain dimensions
     call wrap_inq_dimid  (fid, 'ni' , dimid)
     call wrap_inq_dimlen (fid, dimid, nx   )
     call wrap_inq_dimid  (fid, 'nj' , dimid)
     call wrap_inq_dimlen (fid, dimid, ny   )
     call wrap_inq_dimid  (fid, 'nv' , dimid)
     call wrap_inq_dimlen (fid, dimid, nv   )

     if (nx/=(imt_global-2)) then
        write(6,*)"nx=",nx,"imt=",imt_global
        stop
     end if
  !   if (ny/=(jmt_global-1)) then
  !      write(6,*)"ny=",ny,"jmt==",jmt_global
  !      stop
  !   end if
  end if


  call shr_mpi_bcast(nx,mpi_comm_ocn,"nx")
  call shr_mpi_bcast(ny,mpi_comm_ocn,"ny")
  call shr_mpi_bcast(nv,mpi_comm_ocn,"nv")

     ! obtain grid variables
     allocate(xc(nx,ny))
     allocate(yc(nx,ny))
     allocate(xv(nv,nx,ny))
     allocate(yv(nv,nx,ny))
     allocate(mask(nx,ny))
     allocate(area(nx,ny))

  if (mytid==0 ) then
     call wrap_inq_varid(fid, 'xc' ,vid)
     call wrap_get_var_realx(fid, vid, xc)

     call wrap_inq_varid(fid, 'yc' ,vid)
     call wrap_get_var_realx(fid, vid, yc)

     call wrap_inq_varid(fid, 'xv' ,vid)
     call wrap_get_var_realx(fid, vid, xv)

     call wrap_inq_varid(fid, 'yv' ,vid)
     call wrap_get_var_realx(fid, vid, yv)

     call wrap_inq_varid(fid, 'mask', vid)
     call wrap_get_var_int(fid,vid,mask)

     call wrap_inq_varid(fid, 'area', vid)
     call wrap_get_var_realx(fid,vid,area)

     call wrap_close(fid)
  end if

  ! TODO
  call shr_mpi_bcast(xc,mpi_comm_ocn,"xc")
  call shr_mpi_bcast(yc,mpi_comm_ocn,"yc")
  call shr_mpi_bcast(xv,mpi_comm_ocn,"xv")
  call shr_mpi_bcast(yv,mpi_comm_ocn,"yv")
  call shr_mpi_bcast(mask,mpi_comm_ocn,"mask")
  call shr_mpi_bcast(area,mpi_comm_ocn,"area")
  


  ! lihuimin, 2012.7.17, nx,ny --> imt,jmt
  allocate (t_cpl(imt,jmt))
  allocate (s_cpl(imt,jmt))
  allocate (u_cpl(imt,jmt))
  allocate (v_cpl(imt,jmt))
  allocate (dhdx(imt,jmt))
  allocate (dhdy(imt,jmt))
  allocate (Q   (imt,jmt))

  !allocate (t_cpl(nx,ny))
  !allocate (s_cpl(nx,ny))
  !allocate (u_cpl(nx,ny))
  !allocate (v_cpl(nx,ny))
  !allocate (dhdx(nx,ny))
  !allocate (dhdy(nx,ny))
  !allocate (Q   (nx,ny))

  allocate (taux (imt,jmt))
  allocate (tauy (imt,jmt))
  allocate (netsw(imt,jmt))
  allocate (lat1 (imt,jmt))
  allocate (sen  (imt,jmt))
  allocate (lwup (imt,jmt))
  allocate (lwdn (imt,jmt))
  allocate (melth(imt,jmt))
  allocate (salt (imt,jmt))
  allocate (prec (imt,jmt))
  allocate (evap (imt,jmt))
  allocate (meltw(imt,jmt))
  allocate (roff (imt,jmt))
  allocate (ifrac(imt,jmt)) 
  allocate (patm (imt,jmt)) 

  ! lihuimin, 2012.7.8
  allocate (duu10n(imt,jmt))

  ! lihuimin 2012.6.15, buffs/r used in msg_pass, so comment it
  !allocate (buffs(nx,ny,nsnd))
  !allocate (buffr(nx,ny,nrcv))
!

   ! lihuimin 2012.6.14, commented
!  if (mytid==0) then
!    call msg_pass('init')
!  endif

#endif

!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = 1,JMT
         DO I = 1,IMT
            UB (I,J)= 0.0D0
            VB (I,J)= 0.0D0
            H0 (I,J)= 0.0D0
            UBP (I,J)= 0.0D0
            VBP (I,J)= 0.0D0
            H0P (I,J)= 0.0D0
         END DO
      END DO

 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               U (I,J,K)= 0.0D0
               V (I,J,K)= 0.0D0
               UP (I,J,K)= 0.0D0
               VP (I,J,K)= 0.0D0
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KMP1
         DO J = 1,JMT
            DO I = 1,IMT
               WS (I,J,K)= 0.0D0
            END DO
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = 1,JMT
         DO I = 1,IMT
            H0L (I,J)= 0.0D0
            H0F (I,J)= 0.0D0
            H0BL (I,J)= 0.0D0
            H0BF (I,J)= 0.0D0
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               UTL (I,J,K)= 0.0D0
               UTF (I,J,K)= 0.0D0
               VTL (I,J,K)= 0.0D0
               VTF (I,J,K)= 0.0D0
            END DO
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,NTRA
      DO J = 1,JMT
         DO I = 1,IMT
            NET (I,J,K)= 0.0D0
         END DO
      END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = 1,JMT
         DO I = 1,IMT
            ITICE (I,J)= 0D0
            ALEAD (I,J)= 0.0D0
            TLEAD (I,J)= 0.0D0
            HI (I,J)= 0.0D0
         END DO
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = 1,JMT
         DO I = 1,IMT
            PXB (I,J)= 0.0D0
            PYB (I,J)= 0.0D0
            PAX (I,J)= 0.0D0
            PAY (I,J)= 0.0D0
            WHX (I,J)= 0.0D0
            WHY (I,J)= 0.0D0
            WGP (I,J)= 0.0D0
         END DO
      END DO
 
!     ------------------------------------------------------------------
!     Output Arrays
!     ------------------------------------------------------------------
      CALL YY00
!
 
      MONTH = 1
 
      IF (NSTART == 1) THEN
      number_day = 1
 
!     ------------------------------------------------------------------
!     READ LEVITUS ANNUAL MEAN TEMPERATURE AND SALINITY
!     ------------------------------------------------------------------
#ifdef SPMD
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      if (mytid==0) then
       iret=nf_open('TSinitial',nf_nowrite,ncid)
       call check_err (iret)
      end if

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
     
      do k=1,km !km
!
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret=nf_get_vara_real(ncid,   5,start,count, buffer_real4)
        call check_err (iret)
       end if
!
      call global_distribute_real(buffer_real4,at(1,1,k,1))
!
       if (mytid == 0 ) then
        iret=nf_get_vara_real(ncid,   6,start,count, buffer_real4)
        call check_err (iret)
       end if
        call global_distribute_real(buffer_real4,at(1,1,k,2))
!Yu
       end do !km
#else
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('TSinitial',nf_nowrite,ncid)
      call check_err (iret)

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1

      iret=nf_get_vara_real(ncid,   5,start,count, at_io(1,1,1,1))
      call check_err (iret)

      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1
      iret=nf_get_vara_real(ncid,   6,start,count, at_io(1,1,1,2))
      call check_err (iret)

      iret = nf_close (ncid)
      call check_err (iret)

     at=at_io
#endif
 
!----------------------------------------------------
!   assign 0 to land grids of TSinital
!----------------------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,IMT
                  AT (I,J,K,1) = AT (I,J,K,1)*VIT(I,J,K)
!                  if(vit(i,j,k)<0.5) at(i,j,k,1)=1.0e30 !linpf 2012Jul26 
                  AT (I,J,K,2) = (AT (I,J,K,2)- 35.0D0)*0.001D0*VIT(I,J,K)
!                  if(vit(i,j,k)<0.5) at(i,j,k,2)=1.0e30 
               END DO
            END DO
         END DO
!
 
         DO N = 1,NTRA
!$OMP PARALLEL DO PRIVATE (K,J,I)
            DO K = 1,KM
               DO J = 1,JMT
                  DO I = 1,IMT
                     ATB (I,J,K,N) = AT (I,J,K,N)
#if (defined BOUNDARY)
                     RESTORE (I,J,K,N) = AT (I,J,K,N)
#endif
                  END DO
               END DO
            END DO
 
 
!nick
!$OMP PARALLEL DO PRIVATE (J,I)
            DO J = 1,JMT
               DO I = 1,IMT
                  ATB (I,J,0,N) = 0.0D0
               END DO
            END DO
!nick
         END DO

      ELSE
 
!     ------------------------------------------------------------------
!     READ INTERMEDIATE RESULTS (fort.22/fort.21)
!     ------------------------------------------------------------------
 
#if (defined BOUNDARY)
#ifdef SPMD
         if (mytid==0) then
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('TSinitial',nf_nowrite,ncid)
      call check_err (iret)

          end if
!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
     
      do k=1,km
!
      if (mytid == 0) then
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=k ; count(3)=1
      start(4)=1 ; count(4)=1

      iret=nf_get_vara_real(ncid,   5,start,count, buffer_real4)
      call check_err (iret)
      end if
      call global_distribute_real(buffer_real4,at(1,1,k,1))
!
      if (mytid == 0 ) then
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=k ; count(3)=1
      start(4)=1 ; count(4)=1
      iret=nf_get_vara_real(ncid,   6,start,count, buffer_real4)
      call check_err (iret)
      end if
      call global_distribute_real(buffer_real4,at(1,1,k,2))
!
       end do
!
      if (mytid == 0 ) then
      iret = nf_close (ncid)
      call check_err (iret)
      end if
!!!!!!!!!!!!!!!!!!!!!!!
#else
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('TSinitial',nf_nowrite,ncid)
      call check_err (iret)

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1
            
      iret=nf_get_vara_real(ncid,   5,start,count, restore_io(1,1,1,1))
      call check_err (iret)
      
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1
      iret=nf_get_vara_real(ncid,   6,start,count, restore_io(1,1,1,2))
      call check_err (iret)
!     
      iret = nf_close (ncid)
      call check_err (iret)
      restore=restore_io
#endif
 
!----------------------------------------------------
!   assign 0 to land grids of TSinital
!----------------------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,IMT
                  RESTORE (I,J,K,1) = at (I,J,K,1)*VIT(I,J,K) 
                  RESTORE (I,J,K,2) = (at (I,J,K,2) - 35.0)*0.001*VIT(I,J,K) 
               END DO
            END DO
         END DO
!
#endif
!
         if (mytid==0) then
         open(22,file=trim(out_dir)//fname,form='unformatted')
         end if
#ifdef SPMD
         if (mytid==0) then
         READ (22)buffer
         end if
         call global_distribute(buffer,h0)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
         call global_distribute(buffer,u(1,1,k))
         end do
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
         call global_distribute(buffer,v(1,1,k))
         end do
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
         call global_distribute(buffer,at(1,1,k,1))
         end do
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
         call global_distribute(buffer,at(1,1,k,2))
         end do
!lhl20110728 for ws
         do k=1,km
         if (mytid==0) READ (22)buffer
         end do
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
!lhl20110728
         if (mytid == 0) then
          read(22)number_month,number_day
           month= number_month
         endif
!         if(number_month.ne.mon0.and.number_day.ne.iday)then 
            write(*,*) 'number_month =',number_month,'mon0=',mon0,&
                       'number_day=',number_day,'iday=',iday
!            write(*,*) 'initial month and day error'
!          stop 
!         endif
#ifdef COUP
         if (nstart==2) then
            month=(cdate/10000-1)*12+mod(cdate,10000)/100
!M
!$OMP PARALLEL DO PRIVATE (J,I)            
            do j=1,jmt
            do i=1,imt
               ! lihuimin, TODO, to be considered
               ! lihuimin, 2012.7.23, coordinate with flux_cpl, ft. yu
               !t_cpl (i,j)  = 273.15+at(i,jmt-j+1,1,1)
               !s_cpl (i,j)  = at(i,jmt-j+1,1,2)*1000.+35.
               t_cpl (i,j)  = 273.15+at(i,j,1,1)
               s_cpl (i,j)  = at(i,j,1,2)*1000.+35.
               ! modi end
               q     (i,j)  = 0.0
               u_cpl (i,j)  = 0.0
               v_cpl (i,j)  = 0.0
               dhdx  (i,j)  = 0.0
               dhdy  (i,j)  = 0.0
            end do
            end do
         else
!LPF 20120815
!for t_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
           call global_distribute(buffer,t_cpl)
!for s_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
           call global_distribute(buffer,s_cpl)
!for u_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
           call global_distribute(buffer,u_cpl)
!for v_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
           call global_distribute(buffer,v_cpl)

!for dhdx 
          if (mytid==0) then
           READ (22)buffer
          end if
           call global_distribute(buffer,dhdx)
!for dhdy 
          if (mytid==0) then
           READ (22)buffer
          end if
           call global_distribute(buffer,dhdy)
!for q
          if (mytid==0) then
           READ (22)buffer
          end if
           call global_distribute(buffer,q)
!            read(22)t_cpl,s_cpl,u_cpl,v_cpl,dhdx,dhdy,q
!            cdate  =  10000*((month-1)/12+1)+100*(mod(month-1,12)+1)+1
!LPF 20120815
         end if
#endif
!         end if
!LPF 20120815
!Yu
      call mpi_bcast(month,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(number_month,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(number_day,1,mpi_integer,0,mpi_comm_ocn,ierr)

#else
         READ (22)H0
         do k=1,km
         read(22) ((u(i,j,k),i=1,imt),j=1,jmt)
         end do
         do k=1,km
         read(22) ((v(i,j,k),i=1,imt),j=1,jmt)
         end do
         do k=1,km
         read(22) ((at(i,j,k,1),i=1,imt),j=1,jmt)
         end do
         do k=1,km
         read(22) ((at(i,j,k,2),i=1,imt),j=1,jmt)
         end do
!lhl20110728
         do k=1,km
         read(22) ((ws(i,j,k),i=1,imt),j=1,jmt)
         end do
!lhl20110728
         read(22)number_month,number_day
         month= number_month
#ifdef COUP
         if (nstart==2) then
            month=(cdate/10000-1)*12+mod(cdate,10000)/100
!$OMP PARALLEL DO PRIVATE (J,I) 
            ! lihuimin, 2012.7.23
            do j=1,jmt!ny
            do i=1,imt!nx
               !t_cpl (i,j)  = 273.15+at(i,jmt_global-j+1,1,1)
               !s_cpl (i,j)  = at(i,jmt_global-j+1,1,2)*1000.+35.
               t_cpl (i,j)  = 273.15+at(i,j,1,1)
               s_cpl (i,j)  = at(i,j,1,2)*1000.+35.
               q     (i,j)  = 0.0
               u_cpl (i,j)  = 0.0
               v_cpl (i,j)  = 0.0
               dhdx  (i,j)  = 0.0
               dhdy  (i,j)  = 0.0
            end do
            end do
         else
            read(22)t_cpl,s_cpl,u_cpl,v_cpl,dhdx,dhdy,q
            cdate  =  10000*((month-1)/12+1)+100*(mod(month-1,12)+1)+1
         end if
#endif
#endif
         CLOSE(22)
 
 
         NMFF = MOD (MONTH -1,12)
 
         CALL VINTEG (U,UB)
         CALL VINTEG (V,VB)
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
!YU         DO J = 2,JMM
            DO J = jst,jmt
               DO I = 1,IMT
                  UP (I,J,K) = U (I,J,K)
                  VP (I,J,K) = V (I,J,K)
                  UTF (I,J,K) = U (I,J,K)
                  VTF (I,J,K) = V (I,J,K)
                  ATB (I,J,K,1) = AT (I,J,K,1)
                  ATB (I,J,K,2) = AT (I,J,K,2)
               END DO
            END DO
         END DO
 
 
!$OMP PARALLEL DO PRIVATE (J,I)
!YU      DO J = 2,JMM
         DO J = jst,jmt
            DO I = 1,IMT
               H0P (I,J)= H0 (I,J)
               UBP (I,J)= UB (I,J)
               VBP (I,J)= VB (I,J)
               H0F (I,J)= H0 (I,J)
               H0BF (I,J)= H0 (I,J)
            END DO
         END DO
      END IF
!
      if (mytid==0)then
          write(6,*)"END-----------INIRUN !"
#ifdef COUP
          call shr_sys_flush(6)
#endif
      endif 

!     do k=1,42
!     do j=1,jmt
!     do i=1,imt
!             if ( j_global(j) > 922 .and. j_global(j) < 930 .and. i_global(i) > 1495.and. i_global(i) < 1508 ) then
!                  if ( at(i,j,k,1) < 0.1D0) then
!                     at(i,j,k,1)= restore(i,j,k,1)
!                     at(i,j,k,2)= restore(i,j,k,2)
!                  end if
!              end if
!     end do
!     end do
!     end do



      deallocate(buffer,buffer_real4)

      RETURN

!         open(999,file='fort22.dat',form='unformatted',status='unknown')
!         write(*,*)'Start to Read'
!         write(999,*)H0_io,U_io,V_io,AT_io,HI_io,ITICE_io,ALEAD_io,MONTH
!         close(999)
      END SUBROUTINE INIRUN
 

