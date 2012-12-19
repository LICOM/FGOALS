module cloud_diagnostics

!---------------------------------------------------------------------------------
! Purpose:
!
! Put cloud physical specifications on the history tape
!  Modified from code that computed cloud optics
!  This is a stub so that new radiation package can work
!
! Author: Byron Boville  Sept 06, 2002
!  Modified Oct 15, 2008
!    
!
!---------------------------------------------------------------------------------

   implicit none
   private
   save

   public :: cloud_diagnostics_init, put_cloud_diagnostics

contains

subroutine cloud_diagnostics_init()
    return
end subroutine cloud_diagnostics_init

subroutine put_cloud_diagnostics(state, pbuf)
    use physics_types, only: physics_state
    use phys_buffer,   only: pbuf_size_max, pbuf_fld
    type(physics_state), intent(in)  :: state        ! state variables
    type(pbuf_fld),      intent(in), dimension(pbuf_size_max) :: pbuf
    return
end subroutine put_cloud_diagnostics

end module cloud_diagnostics
