module inital
  !-----------------------------------------------------------------------
  !BOP
  ! !MODULE:  inital --- Define initial conditions for first run of case
  !
  !----------------------------------------------------------------------- 
  ! !USES:

  implicit none

  private   ! By default everything private to this module

  ! !PUBLIC MEMBER FUNCTIONS:

  public cam_initial   ! Cam initialization (formally inital)
  !
  ! !DESCRIPTION: Module for CAM initialization
  !
  !EOP
  !-----------------------------------------------------------------------

contains

  !
  !----------------------------------------------------------------------- 
  !

  !-----------------------------------------------------------------------
  !BOP
  !ROUTINE:  inital --- Define initial conditions for first run of case
  !
  ! !INTERFACE:

  subroutine cam_initial( dyn_in, dyn_out, NLFileName )
    use cam_logfile, only : iulog
    use chem_surfvals,        only : chem_surfvals_init
    use dyn_comp, only : dyn_init1, dyn_init2, dyn_import_t, dyn_export_t
    use startup_initialconds, only : setup_initial, initial_conds
    use phys_grid, only : phys_grid_init
    ! modules from homme
    use parallel_mod, only : par, initmp
    use namelist_mod, only: readnl
    use dp_coupling, only : d_p_coupling, p_d_coupling
    use pio_types, only : file_desc_t
    use spmd_dyn, only : spmd_readnl

    implicit none

    type(dyn_import_t), intent(out)  :: dyn_in
    type(dyn_export_t), intent(out)  :: dyn_out

    character(len=*), intent(in) :: NLFileName
    integer :: npes_homme

    call initcom()
    !
    ! these need to be called before setup_initial
    !
    !
    ! Read in the number of tasks to be assigned to Homme (needed by initmp)
    call spmd_readnl(NLFileName, npes_homme)
    ! Initialize the homme structure that holds the MPI decomposition information
    par=initmp(npes_homme)

    ! Read the homme specific part of the namelist
    call readnl(par, NLFileName)

    call setup_initial()

    call dyn_init1(NLFileName, par, dyn_in, dyn_out)


    !
    ! Define physics data structures
    !
    if(par%masterproc  ) write(iulog,*) 'Running phys_grid_init()'
    call phys_grid_init( )

    !
    ! Initialize ghg surface values before default initial distributions
    ! are set in inidat.
    !
    call chem_surfvals_init()


    if(par%masterproc  ) write(iulog,*) 'Reading initial data'
    call initial_conds(dyn_in)
    if(par%masterproc  ) write(iulog,*) 'Done Reading initial data'

    call dyn_init2(par, dyn_in%elem)



  end subroutine cam_initial
end module inital
