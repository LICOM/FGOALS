!=================================================================================
! utility routines for the offline radiation driver 
! Francis Vitt -- Created 15 Dec 2009
!=================================================================================
module rad_driver_routines

  use shr_kind_mod,    only: r8=>SHR_KIND_R8, cl=>SHR_KIND_CL, cs=>SHR_KIND_CS

  implicit none
  private

  public :: rad_driver_init
  public :: rad_driver_final

  ! private data:

  character(len=cl) :: rad_drv_infile = 'rad_drv_infile'
  character(len=cs) :: rad_drv_case = 'RAD_DRIVER'

contains

  subroutine rad_driver_init(cam_in, cam_out, phys_state, landm)

    use dyn_comp,        only: dyn_import_t, dyn_export_t
    use physics_types,   only: physics_tend
    use filenames,       only: caseid
    use inital,          only: cam_initial
    use ppgrid,          only: pcols
    use ppgrid,          only: begchunk, endchunk
    use phys_buffer,     only: pbuf
    use physconst,       only: epsilo, latvap, latice, rh2o, cpair, tmelt 
    use spmd_utils,      only: spmdinit
    use cam_history,     only: intht, init_masterlinkedlist
    use cam_pio_utils,   only: init_pio_subsystem
    use runtime_opts,    only: read_namelist
    use camsrfexch_types,only: hub2atm_alloc, atm2hub_alloc
    use physics_types,   only: physics_type_alloc, physics_state_set_grid
    use ESMF_Mod,        only: ESMF_Initialize
    use history_defaults,only: bldfld
    use rad_constituents,only: rad_cnst_init
    use aer_rad_props,   only: aer_rad_props_init
    use cloud_rad_props, only: cloud_rad_props_init
    use radiation,       only: radiation_init
    use shr_orb_mod,     only: shr_orb_params
    use cam_control_mod, only: lambm0, obliqr, eccen, mvelpp
    use camsrfexch_types,only: cam_out_t, cam_in_t
    use physics_types,   only: physics_state
    use rad_data_input,  only: init_rad_data_input
    use rad_solar_var,   only: rad_solar_var_init
    use solar_data,      only: solar_data_init
    use seq_io_mod,      only: seq_io_init1,seq_io_init2
    use seq_comm_mct,    only: seq_comm_init 

    implicit none
    
    !args 
    
    type(cam_out_t),     pointer :: cam_out(:)       ! Output from CAM to surface
    type(cam_in_t) ,     pointer :: cam_in(:)        ! Merged input state to CAM
    type(physics_state), pointer :: phys_state(:)
    real(r8),            pointer :: landm(:,:)     ! land fraction ramp

    !local vars
    
    type(physics_tend ), pointer :: phys_tend(:)

    type(dyn_import_t) :: dyn_in   ! Dynamics import container
    type(dyn_export_t) :: dyn_out  ! Dynamics export container

    integer  :: lchnk
    integer  :: orb_yr
    real(r8) :: mvelp, obliq

    integer  :: Global_Comm
    character(*), parameter :: NLFileName = "drv_in"  ! input namelist filename
    integer, pointer :: atm_petlist(:), lnd_petlist(:), ice_petlist(:), ocn_petlist(:), glc_petlist(:)

#include <mpif.h>

    call spmdinit( MPI_COMM_WORLD )

    Global_Comm=MPI_COMM_WORLD
    call seq_io_init1(NLFileName, Global_Comm)

    if (Global_Comm /= MPI_COMM_NULL) then
       call seq_comm_init(Global_Comm, NLFileName, atm_petlist=atm_petlist, lnd_petlist=lnd_petlist, &
            ice_petlist=ice_petlist, ocn_petlist=ocn_petlist, glc_petlist=glc_petlist)
    end if

    call seq_io_init2()

    call init_pio_subsystem('atm_in')
    call ESMF_Initialize()
    !
    ! Read namelist
    !
    call rad_driver_readnl('atm_in')
    caseid = rad_drv_case

    call read_namelist()

    call init_masterlinkedlist()

    call initindx ()

    call init_rad_data_input( rad_drv_infile )

    call init_clock()

    call cam_initial ( dyn_in, dyn_out, NLFileName='atm_in')

    !
    ! Allocate and setup surface exchange data
    !
    call atm2hub_alloc( cam_out )
    call hub2atm_alloc( cam_in )

    call physics_type_alloc(phys_state, phys_tend, begchunk, endchunk)
    do lchnk = begchunk,endchunk
       call physics_state_set_grid(lchnk, phys_state(lchnk))
    end do

    call solar_data_init()
    call rad_solar_var_init()

    ! Initialize rad constituents and their properties
    call rad_cnst_init(pbuf, phys_state)
    call aer_rad_props_init()
    call cloud_rad_props_init()

    call radiation_init()

    call esinti( epsilo, latvap, latice, rh2o, cpair, tmelt )

    call bldfld()

    call intht()

    allocate ( landm( pcols, begchunk:endchunk ))

    orb_yr = 1990
    call shr_orb_params( orb_yr, eccen, obliq, mvelp, obliqr, lambm0, mvelpp, .false. )

  end subroutine rad_driver_init

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

  subroutine rad_driver_final
    use shr_mpi_mod,     only: shr_mpi_finalize
    use rad_data_input,  only: close_rad_data_input
    implicit none
  
    call close_rad_data_input()
    call shr_mpi_finalize('rad_driver_final')

  endsubroutine rad_driver_final

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

  subroutine init_clock()
    use time_manager, only: timemgr_init, tmgr_dtime=>dtime

    use rad_data_input,  only: calendar, dtime, ref_ymd, ref_tod, data_dtime=>dtime, dates, secs, ntimes

    implicit none

    logical :: perpetual_run    ! If in perpetual mode or not
    integer :: perpetual_ymd    ! Perpetual date (YYYYMMDD)
    integer :: start_ymd        ! Start date (YYYYMMDD)
    integer :: start_tod        ! Start time of day (sec)
    integer :: stop_ymd         ! Stop date (YYYYMMDD)
    integer :: stop_tod         ! Stop time of day (sec)

    start_ymd = dates(1)
    start_tod = secs(1)
    stop_ymd = dates(ntimes)
    stop_tod = secs(ntimes)+data_dtime
    perpetual_run =  .false.
    perpetual_ymd =  0

    if (trim(calendar)=='NOLEAP') calendar = 'NO_LEAP'

    call timemgr_init( &
         calendar_in=calendar, &
         start_ymd=start_ymd, &
         start_tod=start_tod, &
         ref_ymd=ref_ymd, &
         ref_tod=ref_tod, &
         stop_ymd=stop_ymd, &
         stop_tod=stop_tod, &
         perpetual_run=perpetual_run, &
         perpetual_ymd=perpetual_ymd  &
         )

     tmgr_dtime = dtime
  endsubroutine init_clock

!--------------------------------------------------------------------------------

  subroutine rad_driver_readnl( nlfile )
    use namelist_utils,only: find_group_name
    use units,         only: getunit, freeunit
    use abortutils,    only: endrun

    implicit none


    ! arguments
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars
    integer :: unitn, ierr

    namelist /rad_drv_nl/ rad_drv_infile, rad_drv_case

    unitn = getunit()
    open( unitn, file=trim(nlfile), status='old' )
    call find_group_name(unitn, 'rad_drv_nl', status=ierr)
    if (ierr == 0) then
       read(unitn, rad_drv_nl, iostat=ierr)
       if (ierr /= 0) then
          call endrun('rad_driver_readnl: ERROR reading namelist')
       end if
    end if
    close(unitn)
    call freeunit(unitn)

  endsubroutine rad_driver_readnl

endmodule rad_driver_routines
