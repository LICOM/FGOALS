module diag_mod
#include <def-undef.h>
!-------------------------------------------------------------------------------
!
! Author: Yongqiang YU  ( 1  Dec, 2003)
!
!-------------------------------------------------------------------------------
use param_mod 
use pconst_mod 
use output_mod 
use tracer_mod
use work_mod 
!
      implicit none
      logical :: diag_msf,diag_bsf,diag_budget,diag_mth
!
      contains
!

!     ===================
      SUBROUTINE msf
!     ===================
!
      RETURN
      END SUBROUTINE MSF

!
!====================================
      SUBROUTINE BAROSF
!====================================
      implicit none
      return
      END SUBROUTINE BAROSF
!

!     =================
      subroutine diag_tracer
!     =================
!     written by liu hai long 2004, jan 
!
      IMPLICIT none
!
!
      RETURN
      END SUBROUTINE diag_tracer
!
!
!========================================
      SUBROUTINE diag_heat_transport
!========================================
!
      implicit none
      RETURN
      END SUBROUTINE diag_heat_transport

end module diag_mod
