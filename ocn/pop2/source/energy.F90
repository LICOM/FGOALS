!  CVS: $Id: energy.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE ENERGY
!     =================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use dyn_mod
use work_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif

#ifdef COUP
use shr_sys_mod
#endif

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0,et0,es0
      REAL (r8) :: EK,EA,EB,ET,ES
      real(r8)::t0,t1,clock0f,nst
      save t0
      integer :: nnn


!---------------------------------------------------------------------
!     TOTAL K.E.
!---------------------------------------------------------------------

      EK = 0.0
      month=(iyfm-1)*12+mon0

!!$OMP PARALLEL DO PRIVATE (K,J,I),reduction(+:EK)
      DO K = 1,KM
         DO J = JSM,JMM
            DO I = 2,IMM
               EK = EK + DZP (K)* DXDYU (J)* VIV (I,J,K) &
                    * (U (I,J,K)* U (I,J,K) + V (I,J,K)* V (I,J,K))
            END DO
         END DO
      END DO

#ifdef SPMD
      call mpi_reduce(ek,ek0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      ek0= ek0*0.5D0* d0
#else
      EK0= EK *0.5D0* D0
#endif

!---------------------------------------------------------------------
!     AVAILABLE P.E.
!---------------------------------------------------------------------

      EA = 0.0
!!$OMP PARALLEL DO PRIVATE (J,I),reduction(+:EA)
      DO J = JSM,JMM
         DO I = 2,IMM
            EA = EA + DXDYT (J)* VIT (I,J,1)* H0 (I,J)* H0 (I,J)
         END DO
      END DO

#ifdef SPMD
      call mpi_reduce(ea,ea0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      ea0= ea0*0.5D0*d0*g
#else
      EA0= EA *0.5D0* D0* G
#endif

!---------------------------------------------------------------------
!     GLOBAL MEAN SST
!---------------------------------------------------------------------

      EB = 0.0
      nst=0
!!$OMP PARALLEL DO PRIVATE (J,I),reduction(+:EB)
      DO J = JSM,JMM
         DO I = 2,IMM
            EB = EB + DYT (J)* SINT (J)* VIT (I,J,1)* AT (I,J,1,1)
         END DO
      END DO

#ifdef SPMD
      call mpi_reduce(eb,eb0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      eb0= eb0/asea
#else
      EB0= EB / ASEA
#endif

!---------------------------------------------------------------------
!     GLOBAL MEAN TEMPERATURE & SALINITY
!---------------------------------------------------------------------

      ET = 0.0D0
      ES = 0.0D0
!!$OMP PARALLEL DO PRIVATE (K,J,I),reduction(+:ET,ES)
      DO K = 1,KM
         DO J = JSM,JMM
            DO I = 2,IMM
               ET = ET + DZP (K)* DXDYT (J)* VIT (I,J,K)* AT(I,J,K,1)
               ES = ES + DZP (K)* DXDYT (J)* VIT (I,J,K)* AT (I,J,K,2)
            END DO
         END DO
      END DO

#ifdef SPMD
      call mpi_reduce(et,et0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(es,es0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      et0= et0/ vsea
      es0= es0/ vsea *1.0D3
#else
      ET0= ET / VSEA
      ES0= ES / VSEA *1.0D3
#endif

!---------------------------------------------------------------------
!     Calculating CPU time for per-day.
!---------------------------------------------------------------------
      if (mytid==0) then
         t1=clock0f()
      endif

      if (mytid==0 )then
      WRITE (6,FMT='(I5,I3,6D25.15)') MONTH,IDAY,EK0,EA0,EB0,ET0,ES0,t1-t0
!      call flush_(6)
      end if

      if (mytid==0) then
         t0=t1
      endif
!
      RETURN
      END SUBROUTINE ENERGY


      SUBROUTINE chk_var3d(var,ch,a)
!     =======================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0
      REAL (r8) :: EK,EA,EB
      real(r8)  :: var(imt,jmt,km)
      integer :: a
      character :: ch*2

      EA = 0.0
      EB = 0.0
      DO k = 1,km
      DO J = jsm,jem
      DO i = 2,imm
         EA=EA+ DZP(K)*DXDYU (J)*VIV(i,j,k)
         EB=EB+ DZP(K)*DXDYT (J)*VIT(i,j,k)
      END DO
      END DO
      END DO


      EK = 0.0
      DO K = 1,km
         DO J = jsm,jem
            DO I = 2,imm
           if(a.eq.1) then
               EK = EK + DZP (K)* DXDYT (J)* VIT (I,J,K)* VAR(I,J,K)*VAR(I,J,K)
           else
               EK = EK + DZP (K)* DXDYU (J)* VIV (I,J,K)* VAR(I,J,K)*VAR(I,J,K)
           endif
            END DO
         END DO
      END DO


#ifdef SPMD
      call mpi_reduce(ek,ek0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(ea,ea0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(eb,eb0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
           if(a.eq.1) then
               EK0 = EK0/EB0
           else
               EK0 = EK0/EA0
           endif
      if (mytid==0 )then
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
      end if
#else
           if(a.eq.1) then
               EK0 = EK/EB
           else
               EK0 = EK/EA
           endif
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
#endif
!
      RETURN
      END SUBROUTINE chk_var3d

      SUBROUTINE chk_var3d1(var,ch,a)
!     =======================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0
      REAL (r8) :: EK,EA,EB
      real(r8)  :: var(imt,jmt,km)
      integer :: a
      character :: ch*2

      EA = 0.0
      EB = 0.0
      DO k = 1,km
      DO J = jsm,jem
      DO i = 2,imm
         EA=EA+ DZP(K)*DXDYU (J)*VIV(i,j,k)
         EB=EB+ DZP(K)*DXDYT (J)*VIT(i,j,k)
      END DO
      END DO
      END DO


      EK = 0.0
      DO K = 1,km
         DO J = jsm,jem
            DO I = 2,imm
           if(a.eq.1) then
               EK = EK + DZP (K)* DXDYT (J)* VIT (I,J,K)* VAR(I,J,K)
           else
               EK = EK + DZP (K)* DXDYU (J)* VIV (I,J,K)* VAR(I,J,K)
           endif
            END DO
         END DO
      END DO


#ifdef SPMD
      call mpi_reduce(ek,ek0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(ea,ea0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(eb,eb0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
           if(a.eq.1) then
               EK0 = EK0/EB0
           else
               EK0 = EK0/EA0
           endif
      if (mytid==0 )then
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
      end if
#else
           if(a.eq.1) then
               EK0 = EK/EB
           else
               EK0 = EK/EA
           endif
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
#endif
!
      RETURN
      END SUBROUTINE chk_var3d1

      SUBROUTINE chk_var1(var,ch,a)
!     =======================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0
      REAL (r8) :: EK,EA,EB
      real(r8)  :: var(imt,jmt,km)
      integer :: a
      character :: ch*2

      EB = 0.0
      DO k = 1,km-1
      DO J = jsm,jem
      DO i = 2,imm
         EA=EA+ DZP(K+1)*DXDYU (J)*VIV(i,j,k+1)
         EB=EB+ DZP(K+1)*DXDYT (J)*VIT(i,j,k+1)
      END DO
      END DO
      END DO

      EK = 0.0
      DO K = 1,km-1
         DO J = jsm,jem
            DO I = 2,imm
           if(a.eq.1) then
               EK = EK + DZP (K+1)* DXDYT (J)* VIT (I,J,K+1)* VAR(I,J,K)
           else
               EK = EK + DZP (K+1)* DXDYU (J)* VIV (I,J,K+1)* VAR(I,J,K)
           endif
            END DO
         END DO
      END DO


#ifdef SPMD
      call mpi_reduce(ek,ek0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(ea,ea0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(eb,eb0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
           if(a.eq.1) then
               EK0 = EK0/EB0
           else
               EK0 = EK0/EA0
           endif
      if (mytid==0 )then
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
      end if
#else
           if(a.eq.1) then
               EK0 = EK/EB
           else
               EK0 = EK/EA
           endif
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
#endif

      RETURN
      END SUBROUTINE chk_var1


      SUBROUTINE chk_var2d(var,ch,a)
!     =======================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif

      IMPLICIT NONE

      REAL (r8) :: EK,EA,EB
      REAL (r8) :: EK0,EA0,EB0
      REAL(r8)  :: var(imt,jmt)
      integer :: a
      character :: ch*2

      EA = 0.0
      EB = 0.0
      DO J = jsm,jem
      DO i = 2,imm
         EA=EA+DYR(J)*SINU(J)*VIV(i,j,1)
         EB=EB+DYT(J)*SINT(J)*VIT(i,j,1)
      END DO
      END DO


      EK = 0.0
         DO J = jsm,jem
            DO I = 2,imm
           if(a.eq.1) then
               EK = EK+DYT(J)*VIT(I,J,1)*VAR(I,J)
           else
               EK = EK+DYR(J)*VIV(I,J,1)*VAR(I,J)
           endif
            END DO
         END DO

#ifdef SPMD
      call mpi_reduce(ek,ek0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(ea,ea0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(eb,eb0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
           if(a.eq.1) then
               EK0 = EK0/EB0
           else
               EK0 = EK0/EA0
           endif
      if (mytid==0)then
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
      end if
#else
           if(a.eq.1) then
               EK0 = EK/EB
           else
               EK0 = EK/EA
           endif
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
#endif

      RETURN
      END SUBROUTINE chk_var2d


