!  CVS: $Id: run_time.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =======================
      SUBROUTINE RUN_TIME(SUBNAME)
!     =======================
 
#include <def-undef.h>
use precision_mod
use param_mod

#ifdef SPMD
use msg_mod
#endif
      IMPLICIT NONE
 
      real(r4) :: t0,t1,clock0f
      character(len=*) :: SUBNAME
      save t0
 
!---------------------------------------------------------------------
!     Calculating CPU time for per-day.
!---------------------------------------------------------------------
#ifdef SPMD
      if (mytid==0) then
#endif
         t1=clock0f()
#ifdef SPMD
      endif
#endif

#ifdef SPMD
      if (mytid==0 )then
      WRITE (6,FMT='(A,E15.7)') SUBNAME,t1-t0
      end if
#else
      WRITE (6,FMT='(A,E15.7)') SUBNAME,t1-t0
#endif

#ifdef SPMD
      if (mytid==0) then
#endif
         t0=t1
#ifdef SPMD
      endif
#endif
      RETURN
      END SUBROUTINE RUN_TIME
 
 
