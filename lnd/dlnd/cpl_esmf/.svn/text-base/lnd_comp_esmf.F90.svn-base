!===============================================================================
! SVN: $Id: lnd_comp_esmf.F90 9228 2008-03-14 00:24:55Z tcraig $ 
! SVN: $URL: https://svn-ccsm-models.cgd.ucar.edu/dlnd7/branch_tags/drva_dlnd7_070824_tags/drva13_dlnd7_071129/lnd_comp_esmf.F90 $
!===============================================================================

module lnd_comp_esmf

   use shr_kind_mod, only:  R8=>SHR_KIND_R8, IN=>SHR_KIND_IN, &
                            CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use shr_sys_mod   ! shared system calls

   use seq_comm_mct, only: LNDID
   use seq_cdata_mod
   use seq_infodata_mod

   use esmfshr_mod

   use dlnd_comp_mod
   use esmf_mod
   use perf_mod
   use mct_mod

   implicit none

   public :: lnd_init_esmf
   public :: lnd_run_esmf
   public :: lnd_final_esmf
   public :: lnd_register_esmf

   private ! except

   type(seq_infodata_type)  :: infodata
   type(seq_cdata)     :: cdata_l
   type(mct_gsMap)     :: gsmap_l
   type(mct_gGrid)     :: ggrid_l
   type(mct_aVect)     :: x2l
   type(mct_aVect)     :: l2x
   type(seq_cdata)     :: cdata_r
   type(mct_gsMap)     :: gsmap_r
   type(mct_gGrid)     :: ggrid_r
   type(mct_aVect)     :: r2x
   type(seq_cdata)     :: cdata_s
   type(mct_gsMap)     :: gsmap_s
   type(mct_gGrid)     :: ggrid_s
   type(mct_aVect)     :: x2s
   type(mct_aVect)     :: s2x

   !----- formats -----
   character(*),parameter :: subName =  "(lnd_comp_esmf) "

   save ! save everything

!
! Author: Fei Liu
! This module is ESMF compliant lnd data component 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================

subroutine lnd_register_esmf(comp, rc)

   implicit none

   type(ESMF_GridComp)  :: comp
   integer, intent(out) :: rc

   rc = ESMF_SUCCESS

   print *, "In lnd register routine"
   ! Register the callback routines.

   call ESMF_GridCompSetEntryPoint(comp, ESMF_SETINIT, lnd_init_esmf, phase=ESMF_SINGLEPHASE, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   call ESMF_GridCompSetEntryPoint(comp, ESMF_SETRUN, lnd_run_esmf, phase=ESMF_SINGLEPHASE, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   call ESMF_GridCompSetEntryPoint(comp, ESMF_SETFINAL, lnd_final_esmf, phase=ESMF_SINGLEPHASE, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

end subroutine

!===============================================================================

subroutine lnd_init_esmf(comp, import_state, export_state, EClock, rc)
   !----------------------------------------------------------
   
   implicit none

   !----- arguments -----
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc
   
   !----- local -----
   character(CL)    :: NLFilename
   type(ESMF_Array) :: Ex2l, El2x, Edoml
   type(ESMF_Array) ::       Er2x, Edomr
   type(ESMF_Array) :: Ex2s, Es2x, Edoms
   integer(IN)      :: phase
   
   character(*),parameter :: subName = "(lnd_init_esmf) "
   character(ESMF_MAXSTR) :: convCIM, purpComp
   !----------------------------------------------------------
   
   rc = ESMF_SUCCESS

   NLFilename = 'unused'

   call esmfshr_infodata_state2infodata(export_state,infodata)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   call seq_infodata_GetData(infodata,lnd_phase=phase)

   if (phase == 1) then
      call seq_cdata_init(cdata_l,LNDID,ggrid_l,gsmap_l,infodata,'dlnd_l')
      call seq_cdata_init(cdata_r,LNDID,ggrid_r,gsmap_r,infodata,'dlnd_r')
      call seq_cdata_init(cdata_s,LNDID,ggrid_s,gsmap_s,infodata,'dlnd_s')
   else
      call ESMF_StateGet(import_state, itemName="x2l", array=Ex2l, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call esmf2mct_copy(Ex2l, x2l, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_StateGet(import_state, itemName="x2s", array=Ex2s, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call esmf2mct_copy(Ex2s, x2s, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   endif

   call dlnd_comp_init(EClock, cdata_l, x2l, l2x, cdata_r, r2x, cdata_s, x2s, s2x, NLFilename)

   call esmfshr_infodata_infodata2state(infodata,export_state,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   if (phase == 1) then
      Edoml = mct2esmf_init(ggrid_l%data,gsmap_l,name='domain_l',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(ggrid_l%data,Edoml,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      Edomr = mct2esmf_init(ggrid_r%data,gsmap_r,name='domain_r',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(ggrid_r%data,Edomr,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      Edoms = mct2esmf_init(ggrid_s%data,gsmap_s,name='domain_s',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(ggrid_s%data,Edoms,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      El2x = mct2esmf_init(l2x,gsmap_l,name='l2x',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(l2x,El2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      Ex2l = mct2esmf_init(x2l,gsmap_l,name='x2l',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      Er2x = mct2esmf_init(r2x,gsmap_r,name='r2x',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(r2x,Er2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      Es2x = mct2esmf_init(s2x,gsmap_s,name='s2x',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(s2x,Es2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      Ex2s = mct2esmf_init(x2s,gsmap_s,name='x2s',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      call ESMF_StateAdd(export_state,Edoml,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_StateAdd(export_state,El2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_StateAdd(import_state,Ex2l,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      call ESMF_StateAdd(export_state,Edomr,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_StateAdd(export_state,Er2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      call ESMF_StateAdd(export_state,Edoms,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_StateAdd(export_state,Es2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_StateAdd(import_state,Ex2s,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   else
      call ESMF_StateGet(export_state, itemName="l2x", array=El2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(l2x,El2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      call ESMF_StateGet(export_state, itemName="r2x", array=Er2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(s2x,Er2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      call ESMF_StateGet(export_state, itemName="s2x", array=Es2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(s2x,Er2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   endif

    convCIM  = "CIM 1.0"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "DLND", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
                           "Climatological Land Data Model", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
                 "The CESM data models perform the basic function of " // &
                 "reading external data, modifying that data, and then " // &
                 "sending it to the driver via standard CESM coupling " // &
                 "interfaces. The driver and other models have no " // &
                 "fundamental knowledge of whether another component " // &
                 "is fully active or just a data model.  In some cases, " // &
                 "data models are prognostic and also receive and use " // &
                 "some data sent by the driver to the data model.  But " // &
                 "in most cases, the data models are not running " // &
                 "prognostically and have no need to receive any data " // &
                 "from the driver.", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Land", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "Name", "Sam Levis", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "EmailAddress", &
                           "slevis@ucar.edu", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
                           convention=convCIM, purpose=purpComp, rc=rc)

   rc = ESMF_SUCCESS

end subroutine lnd_init_esmf

!===============================================================================

subroutine lnd_run_esmf(comp, import_state, export_state, EClock, rc)

   implicit none

   !----- arguments -----
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc

   !----- local -----
   type(ESMF_Array) :: Ex2l, El2x
   type(ESMF_Array) ::       Er2x
   type(ESMF_Array) :: Ex2s, Es2x
   
   character(*),parameter :: subName = "(lnd_run_esmf) "
   !----------------------------------------------------------
   
   rc = ESMF_SUCCESS

   call esmfshr_infodata_state2infodata(export_state,infodata)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   call ESMF_StateGet(import_state, itemName="x2l", array=Ex2l, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   call esmf2mct_copy(Ex2l, x2l, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   call ESMF_StateGet(import_state, itemName="x2s", array=Ex2s, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   call esmf2mct_copy(Ex2s, x2s, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   call dlnd_comp_run(EClock, cdata_l, x2l, l2x, cdata_r, r2x, cdata_s, x2s, s2x)

   call esmfshr_infodata_infodata2state(infodata,export_state,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   call ESMF_StateGet(export_state, itemName="l2x", array=El2x, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   call mct2esmf_copy(l2x,El2x,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   call ESMF_StateGet(export_state, itemName="r2x", array=Er2x, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   call mct2esmf_copy(r2x,Er2x,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   call ESMF_StateGet(export_state, itemName="s2x", array=Es2x, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   call mct2esmf_copy(s2x,Es2x,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   rc = ESMF_SUCCESS

end subroutine lnd_run_esmf

!===============================================================================

subroutine lnd_final_esmf(comp, import_state, export_state, EClock, rc)

   implicit none

   !----- arguments -----
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc

   !----------------------------------------------------------------------------
   ! Finalize routine 
   !----------------------------------------------------------------------------
  
   rc = ESMF_SUCCESS

   call dlnd_comp_final()
  
end subroutine lnd_final_esmf

!===============================================================================

end module lnd_comp_esmf
