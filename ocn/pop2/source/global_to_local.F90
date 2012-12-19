!  CVS: $Id: global_to_local.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      subroutine global_to_local_1d(global,local)
!     ================
!     To transfer global 1-d data to local processor.
 
#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
      IMPLICIT NONE
      real(r8) :: global(jmt_global),local(jmt)

!$OMP PARALLEL DO PRIVATE (J)
      do j=1,jmt
         local(j)=0.0
      end do

!$OMP PARALLEL DO PRIVATE (J)
      do j=1,jmt
         if (j_global(j)<=jmt_global) then
             local(j)=global(j_global(j))
         else
             local(j)=global(jmt_global)
         end if
      end do

#endif

      end subroutine global_to_local_1d
 
!--------------------------------------------------------------
      subroutine global_to_local_2d(global,local)
!     ================
!     To transfer global 2-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
      IMPLICIT NONE
      integer :: j_end
      real(r8) :: global(imt,jmt_global),local(imt,jmt)

      if (mytid==0) then
         do n=1,nproc-1
!
!$OMP PARALLEL DO PRIVATE (J,I)
            do j=1,jmt
            do i=1,imt
               if ((j_start(n+1)+j-1)<=jmt_global) then
                   local(i,j)=global(i,j_start(n+1)+j-1)
               else
                   local(i,j)=global(i,jmt_global)
               end if
            end do
            end do
            call mpi_send(local,imt*jmt,MPI_PR,n,tag_1d,mpi_comm_ocn,ierr)
         end do
!
!$OMP PARALLEL DO PRIVATE (J,I)
         do j=1,jmt
         do i=1,imt
             local(i,j)=global(i,j_global(j))
         end do
         end do
      else
         do n=1,nproc-1
            if (mytid==n) then
               call mpi_recv(local,imt*jmt,MPI_PR,0,tag_1d,mpi_comm_ocn,status,ierr)
            end if
         end do
      end if


#endif

      end subroutine global_to_local_2d

!--------------------------------------------------------------
      subroutine global_to_local_3d(global,local,kk)
!     ================
!     To transfer global 2-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
      IMPLICIT NONE
      integer :: kk,kkk
      real(r8) :: global(imt,jmt_global,kk),local(imt,jmt,kk)

!$OMP PARALLEL DO PRIVATE (KKK,I,J)
      do kkk=1,kk
      do j=1,jmt
      do i=1,imt
         local(i,j,k)=0.0
      end do
      end do
      end do

!$OMP PARALLEL DO PRIVATE (KKK,I,J)
      do kkk=1,kk
      do j=1,jmt
         if (j_global(j)<jmt_global) then
             do i=1,imt
                local(i,j,kkk)=global(i,j_global(j),kkk)
             end do
         else 
             do i=1,imt
                local(i,j,kkk)=global(i,jmt_global,kkk)
             end do
         end if
      end do
      end do

#endif

      end subroutine global_to_local_3d
 
!--------------------------------------------------------------
      subroutine global_to_local_4d(global,local,kk,mm)
!     ================
!     To transfer global 4-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
      IMPLICIT NONE
      integer :: kk,mm,j_end,kkk
      real(r8) :: global(imt,jmt_global,kk,mm),local(imt,jmt,kk,mm)


      if (mytid==0) then
         do n=1,nproc-1
            do m=1,mm
            do kkk=1,kk
!$OMP PARALLEL DO PRIVATE (J,I)
            do j=1,jmt
            do i=1,imt
               if ((j_start(n+1)+j-1)<=jmt_global) then
                   local(i,j,kkk,m)=global(i,j_start(n+1)+j-1,kkk,m)
               else
                   local(i,j,kkk,m)=global(i,jmt_global,kkk,m)
               end if
            end do
            end do
            end do
            end do
            call mpi_send(local,imt*jmt*kk*mm,MPI_PR,n,tag_1d,mpi_comm_ocn,ierr)
         end do
!
         do m=1,mm
         do kkk=1,kk
!$OMP PARALLEL DO PRIVATE (J,I)
         do j=1,jmt
         do i=1,imt
             local(i,j,kkk,m)=global(i,j_global(j),kkk,m)
         end do
         end do
         end do
         end do
      else
         do n=1,nproc-1
            if (mytid==n) then
!
               call mpi_recv(local,imt*jmt*kk*mm,MPI_PR,0,tag_1d,mpi_comm_ocn,status,ierr)
            end if
         end do
      end if

#endif

      end subroutine global_to_local_4d

!--------------------------------------------------------------
      subroutine global_to_local_4d_real(global,local,kk,mm)
!     ================
!     To transfer global 4-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
!chaolee use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nx_proc,ny_proc,status,mpi_comm_ocn
      IMPLICIT NONE
      integer :: kk,mm,i_end,j_end,kkk,n1,n2,nfile
!chaolee      real(r4) :: global(imt,jmt_global,kk,mm)
      real(r4) :: global(imt_global,jmt_global,kk,mm)
      real(r4) :: tmp(imt,jmt,kk,mm)
      real(r8) :: local(imt,jmt,kk,mm)


      if (mytid==0) then
!chaolee         do n=1,nproc-1
       do n1=0,ny_proc-1
       do n2=0,nx_proc-1
            if((n1+n2) .gt. 0) then
            n=n1*nx_proc+n2
!

            do m=1,mm
            do kkk=1,kk
!$OMP PARALLEL DO PRIVATE (J,I)
            do j=1,jmt
            do i=1,imt
               if ((j_start(n+1)+j-1)<=jmt_global) then
                   tmp(i,j,kkk,m)=global(i_start(n+1)+i-1,j_start(n+1)+j-1,kkk,m)
               else
                   tmp(i,j,kkk,m)=global(i_start(n+1)+i-1,jmt_global,kkk,m)
               end if

            end do
            end do
            end do
            end do
            call mpi_send(tmp,imt*jmt*kk*mm,MPI_PR1,n,tag_1d,mpi_comm_ocn,ierr)
            end if
         end do
         end do 
!
         do m=1,mm
         do kkk=1,kk
!$OMP PARALLEL DO PRIVATE (J,I)
         do j=1,jmt
         do i=1,imt
!chaolee             tmp(i,j,kkk,m)=global(i,j_global(j),kkk,m)
             tmp(i,j,kkk,m)=global(i_global(i),j_global(j),kkk,m)
         end do
         end do
         end do
         end do
      else

!chaolee         do n=1,nproc-1
         do n1=0,ny_proc-1
         do n2=0,nx_proc-1
               n=n1*nx_proc+n2
            if (mytid==n .and. n/=0) then
               call mpi_recv(tmp,imt*jmt*kk*mm,MPI_PR1,0,tag_1d,mpi_comm_ocn,status,ierr)
            end if
         end do
         end do
      end if
!
      local=tmp
      call mpi_barrier(mpi_comm_ocn,ierr)
!
#endif

      end subroutine global_to_local_4d_real

!--------------------------------------------------------------
      subroutine global_distribute(global,local)
!     ================
!     To transfer global 4-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
!chaolee use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nx_proc,ny_proc,status,mpi_comm_ocn
      IMPLICIT NONE
!chaolee      integer :: kk,mm,j_start,j_end,kkk
      integer :: i_end,j_end,n1,n2,nfile
!chaolee      real(r4) :: global(imt,jmt_global,kk,mm)
!PY      real(r4) :: global(imt_global,jmt_global)
      real(r8) :: global(imt_global,jmt_global)
!PY      real(r4) :: tmp(imt,jmt)
      real(r8) :: tmp(imt,jmt)
      real(r8) :: local(imt,jmt)


!chaolee         do n=1,nproc-1
!       do n1=0,ny_proc-1
!       do n2=0,nx_proc-1
!             n=n1*nx_proc+n2

      if (mytid==0) then
        do n=1,nx_proc*ny_proc-1
!
!$OMP PARALLEL DO PRIVATE (J,I)
            do j=1,jmt
            do i=1,imt
               if ((j_start(n+1)+j-1)<=jmt_global) then
                   tmp(i,j)=global(i_start(n+1)+i-1,j_start(n+1)+j-1)
               else
                   tmp(i,j)=0.0
               end if
            end do
            end do
            call mpi_send(tmp,imt*jmt,MPI_PR,n,tag_3d,mpi_comm_ocn,ierr)
         end do
!

!$OMP PARALLEL DO PRIVATE (J,I)
         do j=1,jmt
         do i=1,imt
             tmp(i,j)=global(i_global(i),j_global(j))
         end do
         end do

      else

          do n=1,nx_proc*ny_proc-1 
            if (mytid==n .and. n/=0) then
               call mpi_recv(tmp,imt*jmt,MPI_PR,0,tag_3d,mpi_comm_ocn,status,ierr)
            end if
         end do
      end if
!
!     call mpi_barrier(mpi_comm_ocn,ierr)
      local=tmp
!
#endif

      end subroutine global_distribute

      subroutine global_distribute_real(global,local)
!     ================
!     To transfer global 4-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
!chaolee use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nx_proc,ny_proc,status,mpi_comm_ocn
      IMPLICIT NONE
!chaolee      integer :: kk,mm,j_start,j_end,kkk
      integer :: i_end,j_end,n1,n2,nfile
!chaolee      real(r4) :: global(imt,jmt_global,kk,mm)
      real(r4) :: global(imt_global,jmt_global)
      real(r8) :: tmp(imt,jmt)
      real(r8) :: local(imt,jmt)


      if (mytid==0) then
        do n=1,nx_proc*ny_proc-1
            if(n .gt. 0) then
!
!$OMP PARALLEL DO PRIVATE (J,I)
            do j=1,jmt
            do i=1,imt
               if ((j_start(n+1)+j-1)<=jmt_global) then
                   tmp(i,j)=global(i_start(n+1)+i-1,j_start(n+1)+j-1)
               else
                   tmp(i,j)=0.0
               end if
            end do
            end do

            call mpi_send(tmp,imt*jmt,MPI_PR,n,tag_3d,mpi_comm_ocn,ierr)
            end if
         end do
!

!$OMP PARALLEL DO PRIVATE (J,I)
         do j=1,jmt
         do i=1,imt
             tmp(i,j)=global(i_global(i),j_global(j))
         end do
         end do

      else

          do n=1,nx_proc*ny_proc-1 
            if (mytid==n .and. n/=0) then
               call mpi_recv(tmp,imt*jmt,MPI_PR,0,tag_3d,mpi_comm_ocn,status,ierr)
            end if
         end do
      end if
!
      local=tmp
!     call mpi_barrier(mpi_comm_ocn,ierr)
!
#endif

      end subroutine global_distribute_real
!  CVS: $Id: global_to_local.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!
      subroutine distribute_1d(local,global,kk)
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
      real(r8) :: global(imt_global,jmt,kk),local(imt,jmt,kk)
!
      ix=mod(mytid,nx_proc)
      iy=(mytid-ix)/nx_proc
!
      if (ix == 0) then
         do n=mytid+1, mytid+nx_proc-1
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
            do k=1,kk
            do j=1,jmt
            do i=1,imt
                local(i,j,k)=global(i+(n-mytid)*(imt-2),j,k)
            end do
            end do
            end do
            call mpi_send(local,imt*jmt*kk,MPI_PR,n,tag_1d,mpi_comm_ocn,ierr)
         end do
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
         do k=1,kk
         do j=1,jmt
         do i=1,imt
             local(i,j,k)=global(i,j,k)
         end do
         end do
         end do
      else
         do n=iy*nx_proc+1, (iy+1)*nx_proc
            if ( mytid == n) then
               call mpi_recv(local,imt*jmt*kk,MPI_PR,iy*nx_proc,tag_1d,mpi_comm_ocn,status,ierr)
            end if
         end do
      end if
#endif

      end subroutine distribute_1d
