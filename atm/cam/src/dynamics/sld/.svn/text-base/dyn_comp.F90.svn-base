module dyn_comp
!----------------------------------------------------------------------- 
! 
! Purpose: contains calls to setup defualt history stuff that has not found
!          a proper home yet. Shouldn't really exist.
!
! Public functions/subroutines:
!   bldfld
! 
! Author: B.A. Boville from code in cam_history.F90
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use constituents, only: pcnst, cnst_name, cnst_longname
  use constituents, only: tendnam, fixcnam, tottnam, hadvnam, vadvnam
  use ppgrid,       only: pver, pverp
  use pmgrid,       only: plev, plevp, dyndecomp_set

  use cam_history,  only: dyn_decomp, addfld, add_default

  implicit none

  PRIVATE

  public :: dyn_init, dyn_import_t, dyn_export_t

  ! these structures are not used in this dycore, but are included
  ! for source code compatibility.  
  type dyn_import_t
     integer :: placeholder
  end type dyn_import_t
  type dyn_export_t
     integer :: placeholder
  end type dyn_export_t

CONTAINS

!#######################################################################
  subroutine dyn_init(nlfilename)
!
!----------------------------------------------------------------------- 
! 
! Purpose: 
!
! Build Master Field List of all possible fields in a history file.  Each field has 
! associated with it a "long_name" netcdf attribute that describes what the field is, 
! and a "units" attribute.
! 
! Method: Call a subroutine to add each field
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

  use spmd_utils,     only: masterproc
  use sld_control_mod, only : dyn_sld_readnl
#if (defined SPMD)
  use spmd_dyn,        only : spmd_readnl,spmdinit_dyn
#endif

! ARGUMENTS:
    character(len=*), intent(in) :: nlfilename

!
! Local workspace
!
    integer m                     ! Index
    dyndecomp_set = .true.

    call dyn_sld_readnl(nlfilename)

#if (defined SPMD)
    call spmd_readnl(nlfilename)
    call spmdinit_dyn()
#endif 

!
! Call addfld to add each field to the Master Field List.
!

!----------------------------------------------------------------------------
! Dynamics variables which belong in dynamics specific initialization modules
!----------------------------------------------------------------------------

    call addfld ('ETADOT  ','1/s ',plevp,'A','Vertical (eta) velocity',dyn_decomp)
    call addfld ('U&IC    ','m/s ',plev, 'I','Zonal wind'                                    ,dyn_decomp )
    call addfld ('V&IC    ','m/s ',plev, 'I','Meridional wind'                               ,dyn_decomp )
    call add_default ('U&IC       ',0, 'I')
    call add_default ('V&IC       ',0, 'I')

    call addfld ('PS&IC      ','Pa      ',1,    'I','Surface pressure'                              ,dyn_decomp )
    call addfld ('T&IC       ','K       ',plev, 'I','Temperature'                                   ,dyn_decomp )
    call add_default ('PS&IC      ',0, 'I')
    call add_default ('T&IC       ',0, 'I')
    do m = 1,pcnst
       call addfld (trim(cnst_name(m))//'&IC','kg/kg   ',plev, 'I',cnst_longname(m)                 ,dyn_decomp )
    end do

    do m = 1,pcnst
       call add_default(trim(cnst_name(m))//'&IC',0, 'I')
    end do
!
! Constituent tracers
!
    do m=1,pcnst
       call addfld (hadvnam(m), 'kg/kg/s ',pver, 'A',trim(cnst_name(m))//' horizontal advection tendency ',dyn_decomp)
       call addfld (vadvnam(m), 'kg/kg/s ',pver, 'A',trim(cnst_name(m))//' vertical advection tendency ',dyn_decomp)
       call addfld (tendnam(m), 'kg/kg/s ',pver, 'A',trim(cnst_name(m))//' total tendency ',dyn_decomp)
       call addfld (tottnam(m), 'kg/kg/s ',pver, 'A',trim(cnst_name(m))//' horz + vert + fixer tendency ',dyn_decomp)
       call addfld (fixcnam(m), 'kg/kg/s ',pver, 'A',trim(cnst_name(m))//' tendency due to slt fixer',dyn_decomp)
    end do
    call addfld ('DUH     ','K/s     ',plev, 'A','U horizontal diffusive heating',dyn_decomp)
    call addfld ('DVH     ','K/s     ',plev, 'A','V horizontal diffusive heating',dyn_decomp)
    call addfld ('DTH     ','K/s     ',plev, 'A','T horizontal diffusive heating',dyn_decomp)
    call addfld ('ENGYCORR','W/m2    ',plev, 'A','Energy correction for over-all conservation',dyn_decomp)
    call addfld ('TFIX    ','K/s     ',1,    'A','T fixer (T equivalent of Energy correction)',dyn_decomp)

    call add_default ('DTH     ', 1, ' ')

    call addfld ('FU      ','m/s2    ',plev, 'A','Zonal wind forcing term',dyn_decomp)
    call addfld ('FV      ','m/s2    ',plev, 'A','Meridional wind forcing term',dyn_decomp)
    call addfld ('UTEND   ','m/s2    ',plev, 'A','U tendency',dyn_decomp)
    call addfld ('VTEND   ','m/s2    ',plev, 'A','V tendency',dyn_decomp)
    call addfld ('TTEND   ','K/s     ',plev, 'A','T tendency',dyn_decomp)
    call addfld ('LPSTEN  ','Pa/s    ',1,    'A','Surface pressure tendency',dyn_decomp)
    call addfld ('VAT     ','K/s     ',plev, 'A','Vertical advective tendency of T',dyn_decomp)
    call addfld ('KTOOP   ','K/s     ',plev, 'A','(Kappa*T)*(omega/P)',dyn_decomp)
!-----------------------------------------------------------------------
! End of dynamics variables
!-----------------------------------------------------------------------

  end subroutine dyn_init

end module dyn_comp
