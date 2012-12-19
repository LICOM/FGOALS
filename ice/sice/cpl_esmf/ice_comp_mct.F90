module ice_comp_mct

! !USES:
  use shr_kind_mod, only: IN=>shr_kind_IN
  use seq_infodata_mod
  use seq_timemgr_mod
  use mct_mod
  use ESMF_MOD

  use seq_cdata_mod

  use esmfshr_mod
  use ice_comp_esmf
!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: ice_init_mct
  public :: ice_run_mct
  public :: ice_final_mct
  public :: ice_register

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

   type(ESMF_GridComp)     :: ice_comp
   type(ESMF_State)        :: import_state, export_state

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine ice_register(ice_petlist, ccsmComp, localComp)

    integer, pointer                  :: ice_petlist(:)
    type(ESMF_CplComp)                :: ccsmComp
    type(ESMF_GridComp),intent(inout) :: localComp

    integer            :: rc

    ice_comp = ESMF_GridCompCreate(name="ice_comp", petList=ice_petlist, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create ice comp')
    call ESMF_GridCompSetServices(ice_comp, ice_register_esmf, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to register ice comp')
    import_state = ESMF_StateCreate(name="ice import", statetype=ESMF_STATE_IMPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create import ice state')
    export_state = ESMF_StateCreate(name="ice export", statetype=ESMF_STATE_EXPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export ice state')

    call ESMF_AttributeLink(ccsmComp, ice_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to link attribute')

    localComp = ice_comp

end subroutine

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ice_init_mct
!
! !DESCRIPTION:
!     initialize stub ice model
!
! !REVISION HISTORY:
!
! !Authors:
!   Fei Liu
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine ice_init_mct( EClock, cdata, x2d, d2x, NLFilename )

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP

    !--- local variables ---
    integer(IN)                           :: mpicom
    integer(IN)                           :: COMPID
    type(mct_gsMap)             , pointer :: gsMap
    type(mct_gGrid)             , pointer :: dom
    type(seq_infodata_type)     , pointer :: infodata
    integer                               :: rc, urc
!-------------------------------------------------------------------------------

    ! Set cdata pointers
    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsMap, dom=dom, infodata=infodata)

    ! Copy infodata to state

    call esmfshr_infodata_infodata2state(infodata,export_state,rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! call into ESMF init method which calls ESMF run method
    call ESMF_GridCompInitialize(ice_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, terminationflag=ESMF_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! copy export_state to infodata
    call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

end subroutine ice_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ice_run_mct
!
! !DESCRIPTION:
!     run method for stub ice model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine ice_run_mct(EClock, cdata, x2d, d2x)

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> stub
    type(mct_aVect)             ,intent(inout) :: d2x        ! stub   -> driver

!EOP
    type(seq_infodata_type)     , pointer :: infodata
    integer                               :: rc, urc

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, infodata=infodata)

!    ! copy export_state to infodata
!    call esmfshr_infodata_infodata2state(infodata, export_state, rc=rc)
!    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call ESMF_GridCompRun(ice_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, terminationflag=ESMF_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

!    ! copy export_state to infodata
!    call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
!    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

end subroutine ice_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ice_final_mct
!
! !DESCRIPTION:
!     finalize method for stub model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine ice_final_mct()

!EOP

   integer                 :: rc, urc

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    call ESMF_GridCompFinalize(ice_comp, importState=import_state, exportState=export_state, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, terminationflag=ESMF_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    ! destroy component and states
    call ESMF_StateDestroy(import_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_StateDestroy(export_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_GridCompDestroy(ice_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

end subroutine ice_final_mct
!===============================================================================

end module ice_comp_mct
