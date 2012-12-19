#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dimensions_mod
#ifdef CAM
  use constituents, only : qsize_d=>pcnst ! _EXTERNAL
#endif
  implicit none
  private

  integer, parameter, public :: nvar = 3 ! FI # dependent variables 
  integer, parameter, public :: np = NP
  integer, parameter, public :: nv = NP

  integer, parameter, public :: npsq = np*np
  integer, parameter, public :: nvsq = nv*nv

  integer, parameter, public :: nlev=PLEV
  integer, parameter, public :: nlevp=nlev+1
#ifndef CAM
  integer, parameter         :: qsize_d=4  ! dimension of scalars array
#endif
  integer  :: qsize=qsize_d        ! actual numer of scalars to advect
                                   ! some test cases can set to zero
  public :: qsize,qsize_d

  integer, parameter, public :: MaxNeighEdges = 8

  integer, public :: ne
  integer, public :: nelem
  integer, public :: nelemd
  integer, public :: nelemdmax
  integer, public :: nPhysProc                          ! This is the number of physics processors/ per dynamics processor
  integer, public :: nnodes,npart,nmpi_per_node
#ifdef _PRIM
  integer, public :: GlobalUniqueCols
#else
  integer, public :: GlobalUniqueColsP, GlobalUniqueColsV
#endif

end module dimensions_mod

