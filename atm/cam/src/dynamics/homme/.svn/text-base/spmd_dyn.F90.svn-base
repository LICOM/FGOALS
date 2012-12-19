module spmd_dyn

  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: SPMD implementation of CAM homme finite element dynamics.
  ! 
  ! Author: CCM Core Group
  ! Modified: P. Worley, September 2002
  ! 
  !-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use parallel_mod, only: initmp,  par
  use abortutils, only : endrun
  use spmd_utils, only:  masterproc, iam, npes

  implicit none
  private

  integer, allocatable, public :: proc(:)
  public spmd_dyn_defaultopts, spmd_dyn_setopts
  public compute_gsfactors, spmdinit_dyn, spmdbuf, spmd_readnl
  logical,public :: local_dp_map=.true.    ! flag indicates that mapping between dynamics 
  !  and physics decompositions does not require 
  !  interprocess communication
  integer, allocatable, public :: proc_smp_map(:)
  integer, allocatable, public :: nlat_p(:)
  integer,public :: nsmps
  integer, public :: mpicom_homme

  integer, public :: block_buf_nrecs         ! number of local grid points (lon,lat,lev)
                                             !  in dynamics decomposition (including level 0)
  integer, public :: chunk_buf_nrecs         ! number of local grid points (lon,lat,lev)
                                             !  in physics decomposition (including level 0)
                                             ! assigned in phys_grid.F90

CONTAINS
  subroutine spmd_dyn_setopts(npr_yz_in, geopktrans_in,       &
       geopkblocks_in,                                &
       force_2d_in, modcomm_transpose_in,             &
       modcomm_geopk_in, modcomm_gatscat_in,          &
       dyn_alltoall_in, dyn_allgather_in,             &
       dyn_equi_by_col_in,                            &
       dyn_npes_in, dyn_npes_stride_in,               &
       modc_sw_dynrun_in, modc_hs_dynrun_in,          &
       modc_send_dynrun_in, modc_mxreq_dynrun_in,     &
       modc_sw_cdcore_in, modc_hs_cdcore_in,          &
       modc_send_cdcore_in, modc_mxreq_cdcore_in,     &
       modc_sw_gather_in, modc_hs_gather_in,          &
       modc_send_gather_in, modc_mxreq_gather_in,     &
       modc_sw_scatter_in, modc_hs_scatter_in,        &
       modc_send_scatter_in, modc_mxreq_scatter_in,   &
       modc_sw_tracer_in, modc_hs_tracer_in,          &
       modc_send_tracer_in, modc_mxreq_tracer_in,     &
       modc_onetwo_in, modc_tracers_in                )


    !----------------------------------------------------------------------
    ! Purpose: Set runtime options
    ! Author: Art Mirin
    !----------------------------------------------------------------------
    !------------------------------Arguments-------------------------------
    ! yz and xy decompositions (npr_y, npr_z, nprxy_x, nprxy_y)
    integer, intent(in), optional :: npr_yz_in(4)
    ! geopotential method (routines geopk, geopk16, and geopk_d)
    ! 0 for transpose method, 1 for method using semi-global z communication 
    !   with optional 16-byte arithmetic, 2 for method using local
    !   z communication; method 0, method 1 with 16-byte arithmetic and 
    !   method 2 are all bit-for-bit across decompositions; method 0
    !   scales better than method 1 with npr_z, and method 1 is superior 
    !   to method 0 for small npr_z. The optimum speed is attained either 
    !   using method 1 with 8-byte arithmetic (standard for geopk16) or 
    !   method 2 when utilizing the optimal value for the associated 
    !   parameter geopkblocks; see geopk.F90.
    integer, intent(in), optional :: geopktrans_in
    ! number of stages to use in geopotential method geopk_d
    integer, intent(in), optional :: geopkblocks_in
    ! option to force transpose computation for 1D decomp.
    ! the only purpose for invoking this option is debugging
    integer, intent(in), optional :: force_2d_in
    ! mod_comm transpose/geopk/gatscat method
    !   0 for temporary contiguous buffers
    !   1 for mpi derived types
    integer, intent(in), optional :: modcomm_transpose_in, modcomm_geopk_in, &
         modcomm_gatscat_in
    ! Additional mod_comm irregular communication options
    integer, intent(in), optional :: modc_sw_dynrun_in
    logical, intent(in), optional :: modc_hs_dynrun_in
    logical, intent(in), optional :: modc_send_dynrun_in
    integer, intent(in), optional :: modc_mxreq_dynrun_in
    integer, intent(in), optional :: modc_sw_cdcore_in
    logical, intent(in), optional :: modc_hs_cdcore_in
    logical, intent(in), optional :: modc_send_cdcore_in
    integer, intent(in), optional :: modc_mxreq_cdcore_in
    integer, intent(in), optional :: modc_sw_gather_in
    logical, intent(in), optional :: modc_hs_gather_in
    logical, intent(in), optional :: modc_send_gather_in
    integer, intent(in), optional :: modc_mxreq_gather_in
    integer, intent(in), optional :: modc_sw_scatter_in
    logical, intent(in), optional :: modc_hs_scatter_in
    logical, intent(in), optional :: modc_send_scatter_in
    integer, intent(in), optional :: modc_mxreq_scatter_in
    integer, intent(in), optional :: modc_sw_tracer_in
    logical, intent(in), optional :: modc_hs_tracer_in
    logical, intent(in), optional :: modc_send_tracer_in
    integer, intent(in), optional :: modc_mxreq_tracer_in
    integer, intent(in), optional :: modc_onetwo_in
    integer, intent(in), optional :: modc_tracers_in
    ! EUL/SLD-only arguments
    integer, intent(in), optional :: dyn_alltoall_in
    integer, intent(in), optional :: dyn_allgather_in
    logical, intent(in), optional :: dyn_equi_by_col_in
    integer, intent(in), optional :: dyn_npes_in
    integer, intent(in), optional :: dyn_npes_stride_in
    !----------------------------------------------------------------------
    integer color, ierror, ntemp

  end subroutine spmd_dyn_setopts

  subroutine spmd_dyn_defaultopts(npr_yz_out, geopktrans_out,    &
       geopkblocks_out,                                  &
       force_2d_out, modcomm_transpose_out,              &
       modcomm_geopk_out, modcomm_gatscat_out,           &
       dyn_alltoall_out, dyn_allgather_out,              &
       dyn_equi_by_col_out,                              &
       dyn_npes_out, dyn_npes_stride_out,                &
       modc_sw_dynrun_out, modc_hs_dynrun_out,           &
       modc_send_dynrun_out, modc_mxreq_dynrun_out,      &
       modc_sw_cdcore_out, modc_hs_cdcore_out,           &
       modc_send_cdcore_out, modc_mxreq_cdcore_out,      &
       modc_sw_gather_out, modc_hs_gather_out,           &
       modc_send_gather_out, modc_mxreq_gather_out,      &
       modc_sw_scatter_out, modc_hs_scatter_out,         &
       modc_send_scatter_out, modc_mxreq_scatter_out,    &
       modc_sw_tracer_out, modc_hs_tracer_out,           &
       modc_send_tracer_out, modc_mxreq_tracer_out,      &
       modc_onetwo_out, modc_tracers_out                 )
    !----------------------------------------------------------------------
    ! Purpose: Return default runtime options
    ! Author: Art Mirin
    !----------------------------------------------------------------------
    !------------------------------Arguments-------------------------------
    ! yz and xy decompositions
    integer, intent(out), optional :: npr_yz_out(4)
    ! geopotential method (routine geopk, geopk16, or geopk_d)
    integer, intent(out), optional :: geopktrans_out
    ! number of stages to use in geopotential method geopk_d
    integer, intent(out), optional :: geopkblocks_out
    ! option to force transpose computation for 1D decomp.
    integer, intent(out), optional :: force_2d_out
    ! mod_comm transpose method
    integer, intent(out), optional :: modcomm_transpose_out
    ! mod_comm geopk method
    integer, intent(out), optional :: modcomm_geopk_out
    ! mod_comm gather/scatter method
    integer, intent(out), optional :: modcomm_gatscat_out
    ! Additional mod_comm irregular communication options
    integer, intent(out), optional :: modc_sw_dynrun_out
    logical, intent(out), optional :: modc_hs_dynrun_out
    logical, intent(out), optional :: modc_send_dynrun_out
    integer, intent(out), optional :: modc_mxreq_dynrun_out
    integer, intent(out), optional :: modc_sw_cdcore_out
    logical, intent(out), optional :: modc_hs_cdcore_out
    logical, intent(out), optional :: modc_send_cdcore_out
    integer, intent(out), optional :: modc_mxreq_cdcore_out
    integer, intent(out), optional :: modc_sw_gather_out
    logical, intent(out), optional :: modc_hs_gather_out
    logical, intent(out), optional :: modc_send_gather_out
    integer, intent(out), optional :: modc_mxreq_gather_out
    integer, intent(out), optional :: modc_sw_scatter_out
    logical, intent(out), optional :: modc_hs_scatter_out
    logical, intent(out), optional :: modc_send_scatter_out
    integer, intent(out), optional :: modc_mxreq_scatter_out
    integer, intent(out), optional :: modc_sw_tracer_out
    logical, intent(out), optional :: modc_hs_tracer_out
    logical, intent(out), optional :: modc_send_tracer_out
    integer, intent(out), optional :: modc_mxreq_tracer_out
    integer, intent(out), optional :: modc_onetwo_out
    integer, intent(out), optional :: modc_tracers_out
    ! EUL/SLD-only arguments
    integer, intent(out), optional :: dyn_alltoall_out
    integer, intent(out), optional :: dyn_allgather_out
    logical, intent(out), optional :: dyn_equi_by_col_out
    integer, intent(out), optional :: dyn_npes_out
    integer, intent(out), optional :: dyn_npes_stride_out

  end subroutine spmd_dyn_defaultopts

  subroutine spmdinit_dyn ()


    return
  end subroutine spmdinit_dyn

  !========================================================================


  subroutine spmdbuf
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: placeholder for buffer allocation routine 
    ! 
    ! Method: 
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------
    implicit none

    return

  end subroutine spmdbuf
  subroutine compute_gsfactors (numperlat, numtot, numperproc, displs)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Compute arguments for gatherv, scatterv
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: numperlat    ! number of elements per latitude
    !
    ! Output arguments
    !
    integer, intent(out) :: numtot               ! total number of elements (to send or recv)
    integer, intent(out) :: numperproc(0:par%nprocs-1) ! per-PE number of items to receive
    integer, intent(out) :: displs(0:par%nprocs-1)     ! per-PE displacements
    !
    ! Local variables
    !
    integer :: p                    ! index

    !!     numtot = numperlat*numlats

    !!     do p=0,par%nprocs-1
    !!        numperproc(p) = numperlat*nlat_p(p)
    !!     end do

    !!     displs(0) = 0
    !!     do p=1,par%nprocs-1
    !!        displs(p) = displs(p-1) + numperproc(p-1)
    !!     end do

  end subroutine compute_gsfactors

 
  subroutine spmd_readnl(nlfilename, dyn_npes)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use cam_logfile,     only: iulog
    use abortutils,      only: endrun
    use mpishorthand

    implicit none

    character(len=*), intent(in) :: nlfilename
    integer, intent(out) :: dyn_npes
    integer :: ierr           ! error code
    integer :: unitn          ! namelist unit number
    integer :: color, nproc_tmp
    character(len=*), parameter ::  subname = "spmd_readnl"

    logical :: dyn_equi_by_col
    integer :: dyn_alltoall
    integer :: dyn_npes_stride
    integer :: dyn_allgather



!   Note that only dyn_npes is currently used by the homme dycore
    namelist /spmd_dyn_inparm/ dyn_alltoall, &
             dyn_allgather,  &
             dyn_equi_by_col,&
             dyn_npes,       &
             dyn_npes_stride 



    dyn_npes = npes

    if (masterproc) then
       write(iulog,*) 'Read in spmd_dyn_inparm namelist from: ', trim(nlfilename)
       unitn = getunit()
       open( unitn, file=trim(nlfilename), status='old' )
       ! Look for spmd_dyn_inparm group name in the input file.  If found, leave the
       ! file positioned at that namelist group.
       call find_group_name(unitn, 'spmd_dyn_inparm', status=ierr)
       if (ierr == 0) then  ! found spmd_dyn_inparm
          read(unitn, spmd_dyn_inparm, iostat=ierr)  ! read the spmd_dyn_inparm namelist group
          if (ierr /= 0) then
             call endrun( subname//':: namelist read returns an'// &
                  ' error condition for spmd_dyn_inparm' )
          end if
       end if
       if (dyn_npes .lt. 1 .or. dyn_npes .gt. npes) then
          call endrun( subname//':: namelist read returns a'// &
               ' bad value for dyn_npes' )
       endif
       write (iulog,*) 'Homme will use ', dyn_npes, '  tasks'
       close( unitn )
       call freeunit( unitn )
    endif

    call mpibcast (dyn_npes,1,mpiint,0,mpicom)

  end subroutine spmd_readnl



end module spmd_dyn
