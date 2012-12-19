module precision_mod
#include <def-undef.h>
!
#ifdef SPMD
#include <mpif.h>
#endif
!-------------------------------------------------------------------------------
#ifdef D_PRECISION
   integer, parameter, public ::               &
      char_len       = 256                    ,&
      char_len_long  = 512                    ,&
      log_kind       = kind(.true.)           ,&
      int_kind       = kind(1)                ,&
      i4             = selected_int_kind(6)   ,&
      i8             = selected_int_kind(13)  ,&
      r4             = selected_real_kind(6)  ,&
      r6             = selected_real_kind(9) , &
      r8             = selected_real_kind(12) ,&
      r16            = selected_real_kind(24)

#ifdef SPMD
integer*4, parameter :: MPI_PR = MPI_DOUBLE_PRECISION
integer*4, parameter :: MPI_PR1 = MPI_REAL
#endif

#else
   integer, parameter, public ::               &
      char_len       = 256                    ,&
      char_len_long  = 512                    ,&
      log_kind       = kind(.true.)           ,&
      int_kind       = kind(1)                ,&
      i4             = selected_int_kind(6)   ,&
      i8             = selected_int_kind(13)  ,&
      r4             = selected_real_kind(6)  ,&
      r6             = selected_real_kind(9) , &
      r8             = selected_real_kind(12) ,&
      r16            = selected_real_kind(24)

#ifdef SPMD
integer*4, parameter :: MPI_PR = MPI_REAL
integer*4, parameter :: MPI_PR1 = MPI_REAL
#endif

#endif
end module precision_mod

