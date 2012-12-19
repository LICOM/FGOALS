!  CVS: $Id: boundary.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine boundary
!     ================
!     To compute bounary of subdomain for the each processor.
 
#include <def-undef.h>
#ifdef SPMD
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn, &
                   nbrs,my_rank,UPUP,DOWN,RIGHT,LEFT,dims,coords,cartcomm,reorder,periods
use dyn_mod, only : buffer
      IMPLICIT NONE


!
      allocate ( buffer(imt_global,imt_global), i_start(nproc), j_start(nproc) )
!
      if ((imt_global-2) /= (imt-2)*nx_proc) then
         write(6,*) "Error in distributing MPI tasks!"
         stop
      end if
      ix=mod(mytid,nx_proc)
      iy=(mytid-ix)/nx_proc
!
      do i=1,imt
         i_global(i)=ix*(imt-num_overlap)+i
      end do
!
      do j=jst,jmt
         if (iy==0) then
            j_global(j)=jst_global+j-1
         else
            j_global(j)=jst_global-1+iy*(jmt-num_overlap)+j
         end if
      end do
!
      if (mytid == 0) then
          i_start(1)= i_global(1)
          j_start(1)= j_global(1)
          do n=1,nproc-1
          call mpi_recv(j_start(n+1),1,mpi_integer,n,tag_1d,mpi_comm_ocn,status,ierr)
          call mpi_recv(i_start(n+1),1,mpi_integer,n,tag_2d,mpi_comm_ocn,status,ierr)
          end do
       else 
          j_start(mytid+1) =j_global(1)
          i_start(mytid+1) =i_global(1)
          call mpi_send(j_start(mytid+1),1,mpi_integer,0,tag_1d,mpi_comm_ocn,ierr)
          call mpi_send(i_start(mytid+1),1,mpi_integer,0,tag_2d,mpi_comm_ocn,ierr)
       end if
!
     
      if (iy==(ny_proc-1)) then
         if (j_global(jmt)<jmt_global.or.j_global(1)>jmt_global) then
           write(6,*) "ERROR in boundary! "
           write(6,*) "j_global(1),j_global(jmt),jmt_global=",j_global(1),j_global(jmt),jmt_global
           stop
         end if
      end if

      do k = 1,km
      if (mytid .eq. 0) then	
      do j = 1,jmm_global
         do i = 2,imt_global
            buffer(i,j)= vit_global(i-1,j,k)*vit_global(i-1,j+1,k)* &
                         vit_global(i,j,k)* vit_global(i,j+1,k)
         end do
            buffer(1,j)= buffer (imt_global-1,j)
      end do
      do i = 1,imt_global
            buffer(i,jmt_global)= 0.0
      end do
	
      end if
      call global_distribute(buffer,viv(1,1,k))
      end do
!
      deallocate (buffer)

      dims(1) = ny_proc
      dims(2) = nx_proc
      periods(1)= 0
      periods(2)= 1
      call MPI_CART_CREATE(MPI_COMM_OCN, 2, dims, periods, reorder,cartcomm, ierr)
!     call MPI_COMM_RANK(cartcomm, my_rank, ierr)
      call MPI_CART_COORDS(cartcomm, mytid, 2, coords, ierr)
      call MPI_CART_SHIFT(cartcomm, 0, 1, nbrs(UPUP), nbrs(DOWN), ierr)
      call MPI_CART_SHIFT(cartcomm, 1, 1, nbrs(LEFT), nbrs(RIGHT),ierr)


#endif
      end subroutine boundary
 
 
