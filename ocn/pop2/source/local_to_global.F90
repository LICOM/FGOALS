!  CVS: $Id: local_to_global.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      subroutine local_to_global(hist_output,rest_output)
!     ================
!     To transfer global 1-d data to local processor.
 
#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use output_mod
use dyn_mod
use tracer_mod
use work_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn

      IMPLICIT NONE
      integer :: itag
      logical :: hist_output,rest_output 

     if (hist_output) then
!        call local_to_global_4d(z0mon,z0mon_io,1,1)
!        call local_to_global_4d(himon,himon_io,1,1)
!        call local_to_global_4d(hdmon,hdmon_io,1,1)
!        call local_to_global_4d(netmon,netmon_io,2,1)
!        call local_to_global_4d(tsmon,tsmon_io,km,1)
!        call local_to_global_4d(ssmon,ssmon_io,km,1)
!        call local_to_global_4d(usmon,usmon_io,km,1)
!        call local_to_global_4d(vsmon,vsmon_io,km,1)
!        call local_to_global_4d(wsmon,wsmon_io,km,1)
!        call local_to_global_4d(icmon,icmon_io,2,1)
#if (defined SMAG_OUT)
!        call local_to_global_4d(am3mon,am3mon_io,km,1)
#endif
!
!        call local_to_global_4d(axmon,axmon_io,km,2)
!        call local_to_global_4d(aymon,aymon_io,km,2)
!        call local_to_global_4d(azmon,azmon_io,km,2)
!        call local_to_global_4d(dxmon,dxmon_io,km,2)
!        call local_to_global_4d(dymon,dymon_io,km,2)
!        call local_to_global_4d(dzmon,dzmon_io,km,2)
!        call local_to_global_4d(ddymon,ddymon_io,km,2)
#ifdef ISO
!        call local_to_global_4d(axmon_iso,axmon_iso_io,km,2)
!        call local_to_global_4d(aymon_iso,aymon_iso_io,km,2)
!        call local_to_global_4d(azmon_iso,azmon_iso_io,km,2)
!        call local_to_global_4d(dxmon_iso,dxmon_iso_io,km,2)
!        call local_to_global_4d(dymon_iso,dymon_iso_io,km,2)
!        call local_to_global_4d(dzmon_iso,dzmon_iso_io,km,2)
!
!        call local_to_global_4d(aaymon_iso,aaymon_iso_io,km,2)
!        call local_to_global_4d(ddymon_iso,ddymon_iso_io,km,2)
#endif
!        call local_to_global_4d(trendmon,trendmon_io,km,2)
!        call local_to_global_4d(penmon,penmon_io,km,1)
!        call local_to_global_4d(mldmon,mldmon_io,1,1)
!        call local_to_global_4d(akmmon,akmmon_io,km,1)
!        call local_to_global_4d(aktmon,aktmon_io,km,1)
!        call local_to_global_4d(aksmon,aksmon_io,km,1)
!
     end if
!
     if (rest_output) then
!        call local_to_global_4d_double(h0,h0_io,1,1)
!        call local_to_global_4d_double(u,u_io,km,1)
!        call local_to_global_4d_double(v,v_io,km,1)
!        call local_to_global_4d_double(at,at_io,km,2)
     end if

!           end do
#endif
      return
      end subroutine local_to_global
! 
!     ================
     subroutine local_to_global_4d(local,global,kk,mm,toru)
!     ================
!    To transfer global 1-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use output_mod, only: spval
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn

      IMPLICIT NONE
      integer :: itag,kk,mm,toru
      real(r4)    :: local(imt,jmt,kk,mm),global(imt_global,jmt_global,kk,mm), &
                   work(imt,jmt,kk,mm)
!
        if (toru.eq.1) then  !t grid
           do m=1,mm
           do k=1,kk
            do i=1,imt
             do j=1,jmt
              if(vit(i,j,k)<0.5) then 
               local(i,j,k,m)=spval
              else
               local(i,j,k,m)=local(i,j,k,m)/(nmonth (mon0))
              endif
             enddo
            enddo
           end do
           end do
        else  !u grid
           do m=1,mm
           do k=1,kk
            do i=1,imt
             do j=1,jmt
              if(viv(i,j,k)<0.5) then 
               local(i,j,k,m)=spval
              else
               local(i,j,k,m)=local(i,j,k,m)/(nmonth (mon0))
              endif
             enddo
            enddo
           end do
           end do
        endif
!
      if (mytid == 0) then
 
!$OMP PARALLEL DO PRIVATE (I,J,M,k)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              global(i,j_global(j),k,m)=local(i,j,k,m)
           end do
           end do
           end do
           end do

           do m=1,mm
           do k=1,kk
            do i=1,imt_global
             do j=1,jst_global
               global(i,j,k,m)=spval
             enddo
            enddo
           end do
           end do
!
        do n=1,nproc-1
!
           call mpi_recv(work,imt*jmt*kk*mm,mpi_real,n,tag_3d,mpi_comm_ocn,status,ierr)
!$OMP PARALLEL DO PRIVATE (I,J)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              if ((j_start(n+1)+j-1)<=jmt_global) then
                 global(i+i_start(n+1)-1,j+j_start(n+1)-1,k,m)=work(i,j,k,m)
              end if
           end do
           end do
           end do
           end do
         end do
      else
         call mpi_send(local,imt*jmt*kk*mm,mpi_real,0,tag_3d,mpi_comm_ocn,ierr)
      end if
!
      call mpi_barrier(mpi_comm_ocn,ierr)
!
#endif
     return
     end subroutine local_to_global_4d

!     ================
     subroutine local_to_global_4d_double(local,global,kk,mm)
!     ================
!    To transfer global 1-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn

      IMPLICIT NONE
      integer :: itag,kk,mm
      real(r8):: local(imt,jmt,kk,mm),work(imt,jmt,kk,mm)
      real(r8):: global(imt_global,jmt_global,kk,mm)

!
      if (mytid == 0) then
 
!$OMP PARALLEL DO PRIVATE (I,J,M,k)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              global(i,j_global(j),k,m)=local(i,j,k,m)
           end do
           end do
           end do
           end do
!
        do n=1,nproc-1
!
           call mpi_recv(work,imt*jmt*kk*mm,MPI_PR,n,tag_3d,mpi_comm_ocn,status,ierr)
!$OMP PARALLEL DO PRIVATE (I,J)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              if ((j_start(n+1)+j-1)<=jmt_global) then
                 global(i+i_start(n+1)-1,j+j_start(n+1)-1,k,m)=work(i,j,k,m)
              end if
           end do
           end do
           end do
           end do
         end do
      else
         call mpi_send(local,imt*jmt*kk*mm,MPI_PR,0,tag_3d,mpi_comm_ocn,ierr)
      end if
!
      call mpi_barrier(mpi_comm_ocn,ierr)
#endif
     return
     end subroutine local_to_global_4d_double

!     ================
     subroutine global_gather(local,global)
!     ================
!    To transfer global 2-d data from local processors to master processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nx_proc,ny_proc,status,mpi_comm_ocn

      IMPLICIT NONE
      integer :: itag,n1,n2,nfile
      real(r4)    :: local(imt,jmt),global(imt_global,jmt_global),work(imt,jmt)
!
      if (mytid == 0) then
 
!$OMP PARALLEL DO PRIVATE (I,J)

           do j=1,jmt
           do i=1,imt
              global(i_global(i),j_global(j))=local(i,j)
           end do
           end do

!
        do n1=0,ny_proc-1
        do n2=0,nx_proc-1
           n=n1*nx_proc+n2
!
           if(n>0) then
           call mpi_recv(work,imt*jmt,mpi_real,n,itag,mpi_comm_ocn,status,ierr)
!$OMP PARALLEL DO PRIVATE (I,J)

           do j=1,jmt
           do i=1,imt
              if ((j_start(n+1)+j)-1 <=jmt_global) then
                 global(i+i_start(n+1)-1,j+j_start(n+1)-1)=work(i,j)
              end if
           end do
           end do
           end if

         end do
         end do
      else
             call mpi_send(local,imt*jmt,mpi_real,0,itag,mpi_comm_ocn,ierr)
      end if
!
      call mpi_barrier(mpi_comm_ocn,ierr)
!
#endif
     return
     end subroutine global_gather
!  CVS: $Id: local_to_global.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      subroutine gather_1d(local,global,kk)
!     ================
!     To transfer global 2-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
      IMPLICIT NONE
      integer :: j_end,kk
      real(r8) :: global(imt_global,jmt,kk),local(imt,jmt,kk),tmp(imt,jmt,kk)
!
      ix=mod(mytid,nx_proc)
      iy=(mytid-ix)/nx_proc
!
      if (ix == 0) then
            do k=1,kk
            do j=1,jmt
            do i=1,imt
                global(i,j,k) = local(i,j,k)
            end do
                global(imt_global,j,k) = local(2,j,k)
            end do
            end do
         do n=mytid+1, mytid+nx_proc-1
            call mpi_recv(tmp,imt*jmt*kk,MPI_PR,n,tag_1d,mpi_comm_ocn,status,ierr)
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
            do k=1,kk
            do j=1,jmt
            do i=2,imm
                global(i+(n-mytid)*(imt-2),j,k) = tmp(i,j,k)
            end do
            end do
            end do
         end do
!
      else
         do n=iy*nx_proc+1, (iy+1)*nx_proc-1
            if ( mytid == n) then
                call mpi_send(local,imt*jmt*kk,MPI_PR,iy*nx_proc,tag_1d,mpi_comm_ocn,ierr)
            end if 
         end do
      end if
#endif

      end subroutine gather_1d
