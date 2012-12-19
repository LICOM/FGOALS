!  CVS: $Id: msg_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module msg_mod
#include <def-undef.h>
#if ( defined SPMD ) || (defined COUP)
#include <mpif.h>
!-------------------------------------------------------------------------------
!
! Author: Yongqiang YU  ( 12 Nov, 2002)
!
!-------------------------------------------------------------------------------
      integer, parameter :: tag_1d=10, tag_2d=100, tag_3d=200,tag_4d=300
      integer:: nproc
      integer :: status (MPI_STATUS_SIZE) ! Status of message
      integer ,save            :: mpi_comm_ocn
      integer,parameter :: UPUP=1, DOWN=2, LEFT=3, RIGHT=4
      integer my_rank, source, dest, nbrs(4), dims(2),coords(2),cartcomm, periods(2), reorder
#endif
end module msg_mod
