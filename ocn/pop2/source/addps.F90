!  CVS: $Id: addps.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE ADDPS
!     ================
!     COMPENSATING THE LOSS OF GROSS MASS
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
      IMPLICIT NONE
      REAL(r8)    :: ERROR,DH00,error0
      ERROR = 0.0D0
 
!!!$OMP PARALLEL DO PRIVATE (J,I),reduction(+:ERROR)
      DO J = JSM,JMM
         DO I = 2,IMM
            ERROR = ERROR + DYT (J)* SINT (J)* H0 (I,J)* VIT (I,J,1)
         END DO
      END DO
#ifdef SPMD
       call mpi_reduce(error,error0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
       call mpi_bcast(error0,1,MPI_PR,0,mpi_comm_ocn,ierr)
      dh00 = - error0/ asea
#else
      DH00 = - ERROR / ASEA
#endif
 
!$OMP PARALLEL DO PRIVATE (J,I), shared(dh00)
      DO J = JST,JMT
         DO I = 1,IMT
            H0 (I,J)= (H0 (I,J) + DH00)* VIT (I,J,1)
         END DO
      END DO
 
      RETURN
      END SUBROUTINE ADDPS
 
