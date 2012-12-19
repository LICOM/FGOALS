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
! !REVISION HISTORY:
!   05.08.11   Kluzek     Creation
!   05.11.10   Sawyer     Now using dyn_import/export_t containers
!   06.04.13   Sawyer     Removed dependency on prognostics
!
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

! !USES:
   use shr_kind_mod,         only : r8 => shr_kind_r8
   use dyn_comp,             only : dyn_import_t, dyn_export_t, dyn_init
   use phys_grid,            only : phys_grid_init
   use chem_surfvals,        only : chem_surfvals_init
   use startup_initialconds, only : setup_initial, initial_conds
   use dynamics_vars,        only : T_FVDYCORE_STATE   
   use dyn_internal_state,   only : get_dyn_state
!-----------------------------------------------------------------------
!
! Arguments
!
   type(dyn_import_t), intent(out)  :: dyn_in
   type(dyn_export_t), intent(out)  :: dyn_out
   character(len=*), intent(in) :: NLFileName

! !DESCRIPTION:
!
!   Define initial conditions for first run of case
!
! !REVISION HISTORY:
!
!   92.06.01      Bath          Creation from CCM1
!   96.03.01      Acker         Modifications
!   96.04.01      Boville       Reviewed
!   01.06.17      Sawyer        Added call to dynamics_init
!   01.07.12      Sawyer        Added arguments to dynamics_init
!   05.11.10      Sawyer        Added dyn_in, dyn_out, dyn_state to init
!   06.04.13      Sawyer        Removed call to initial_prognostics
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   type (T_FVDYCORE_STATE), pointer :: dyn_state
!
!-----------------------------------------------------------------------
!
   call setup_initial()
   !
   ! Initialize dynamics
   !
   dyn_state => get_dyn_state()

   call dyn_init(dyn_state, dyn_in, dyn_out, NLFileName )
   !
   ! Initialize dynamics grid
   !
   call initcom
   !
   ! Define physics data structures
   !
   call phys_grid_init
   !
   ! Initialize ghg surface values before default initial distributions
   ! are set in inidat.
   !
   call chem_surfvals_init()
   call initial_conds( dyn_in )
!EOC
end subroutine cam_initial

!
!----------------------------------------------------------------------- 
!

end module inital
