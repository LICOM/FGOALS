!================================================================================================
! This is the 'none' chemistry module.
! Most of the routines return without doing anything.
!================================================================================================

module chemistry
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use physics_types,       only: physics_state, physics_ptend
  use ppgrid,              only: begchunk, endchunk, pcols
  use phys_buffer,         only: pbuf_size_max, pbuf_fld

  implicit none
  private
  save
  !
  ! Public interfaces
  !
  public chem_is                        ! identify which chemistry is being used
  public chem_register                  ! register consituents
  public chem_is_active                 ! returns true if this package is active (ghg_chem=.true.)
  public chem_implements_cnst           ! returns true if consituent is implemented by this package
  public chem_init_cnst                 ! initialize mixing ratios if not read from initial file
  public chem_init                      ! initialize (history) variables
  public chem_timestep_init             ! time interpolate chemical loss frequencies
  public chem_timestep_tend             ! interface to tendency computation
  public chem_final
  public chem_write_restart
  public chem_read_restart
  public chem_init_restart
  public chem_readnl                    ! read chem namelist 

  interface chem_write_restart
     module procedure chem_write_restart_bin
     module procedure chem_write_restart_pio
  end interface
  interface chem_read_restart
     module procedure chem_read_restart_bin
     module procedure chem_read_restart_pio
  end interface

  ! Private data

!================================================================================================
contains
!================================================================================================

  logical function chem_is (name)

    character(len=*), intent(in) :: name

    chem_is = .false.
    if (name == 'none' ) then
       chem_is = .true.
    end if

  end function chem_is

!================================================================================================

  subroutine chem_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents for parameterized greenhouse gas chemistry
    ! 
    !-----------------------------------------------------------------------

  end subroutine chem_register

!================================================================================================

  subroutine chem_readnl(nlfile)

    ! args

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input


  end subroutine chem_readnl

!================================================================================================

  function chem_is_active()
    !-----------------------------------------------------------------------
    logical :: chem_is_active
    !-----------------------------------------------------------------------
    chem_is_active = .false.
  end function chem_is_active

!================================================================================================

  function chem_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this package
    ! 
    ! Author: B. Eaton
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: chem_implements_cnst        ! return value

    chem_implements_cnst = .false.

  end function chem_implements_cnst

!===============================================================================

  subroutine chem_init(phys_state)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize parameterized greenhouse gas chemistry
    !          (declare history variables)
    ! 
    !-----------------------------------------------------------------------
    use cam_history,    only: addfld, add_default, phys_decomp
    type(physics_state), intent(in):: phys_state(begchunk:endchunk)

  end subroutine chem_init

!===============================================================================

  subroutine chem_timestep_init(phys_state)
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
         is_perpetual
    type(physics_state), intent(in):: phys_state(begchunk:endchunk)                 


  end subroutine chem_timestep_init

!===============================================================================

  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dt, pbuf, fh2o, fsds )
    use cam_history,      only: outfld
    use camsrfexch_types, only: cam_in_t, cam_out_t
    !-----------------------------------------------------------------------
    !
    ! Arguments:
    !
    real(r8),            intent(in)    :: dt              ! time step
    type(physics_state), intent(in)    :: state           ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend           ! indivdual parameterization tendencies
    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(in)    :: cam_out
    real(r8),            intent(out)   :: fh2o(pcols)     ! h2o flux to balance source from chemistry
    type(pbuf_fld),      intent(in)    :: pbuf(pbuf_size_max)
    real(r8),            intent(in)    :: fsds(pcols)     ! longwave down at sfc

    return
  end subroutine chem_timestep_tend

!===============================================================================

  subroutine chem_init_cnst(name, q, gcid)

    character(len=*), intent(in) :: name         ! constituent name
    real(r8), intent(out) :: q(:,:)   !  mass mixing ratio (gcol, plev)
    integer, intent(in) :: gcid(:)    !  global column id

    return
  end subroutine chem_init_cnst

!===============================================================================
  subroutine chem_final
    return
  end subroutine chem_final
!===============================================================================
  subroutine chem_write_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    return
  end subroutine chem_write_restart_bin
!===============================================================================
  subroutine chem_read_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    return
  end subroutine chem_read_restart_bin
!===============================================================================
  subroutine chem_write_restart_pio( File )
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_write_restart_pio
!===============================================================================
  subroutine chem_read_restart_pio( File )
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_read_restart_pio
!===============================================================================
  subroutine chem_init_restart(File)
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_init_restart

end module chemistry
