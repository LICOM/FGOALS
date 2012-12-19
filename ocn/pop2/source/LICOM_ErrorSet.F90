!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 subroutine LICOM_ErrorSet(errorCode, errorMsg)

! !DESCRIPTION:
!  This routine sets an error code to POP\_Fail and adds a message to
!  the error log for later printing.
!
! !REVISION HISTORY:
!  same as module

   use precision_mod
   use param_mod
   use msg_mod
   use POP_KindsMod

   integer (POP_i4), intent(out) :: &
      errorCode              ! Error code to set to fail

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      errorMsg               ! message to add to error log for printing

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  Set error code to fail
!
!-----------------------------------------------------------------------


   if (mytid == 0) then
       write(6,*) errorMsg
       call flush(6)
   end if
!-----------------------------------------------------------------------
!EOC

 end subroutine LICOM_ErrorSet

