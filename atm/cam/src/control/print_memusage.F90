subroutine print_memusage (string)
!
! Purpose: interface function to C routine print_memusage
!
! Method: Handle all character strings locally to avoid machine-dependent
!         inter-language communication of character variables. 
!         print_memusage will return variables in native types.
!
   use spmd_utils,  only: iam, masterproc
   use cam_logfile, only: iulog
!
! Arguments
!
   character(len=*), intent(in) :: string
!
! Local workspace
!
   integer :: size         ! process size
   integer :: rss          ! process resident set size
   integer :: share        ! process shared memory size
   integer :: text         ! process text size
   integer :: datastack    ! data + stack memory
   integer :: ret          ! return code from get_memusage
!
! Externals
!
   integer, external :: get_memusage

#if defined(BGL) || defined(BGP)
   ret = -1
#else
   ret = get_memusage (size, rss, share, text, datastack)
#endif
   if (masterproc) then
      write(iulog,'(a,i3," ",a,a)')'print_memusage iam ', iam, string, '. -1 in the next line means unavailable'
   endif
   if (ret == 0) then
#ifndef DEBUG
      if (masterproc) &
#endif
      write(iulog,'(a,5i8)')'print_memusage: size, rss, share, text, datastack=', &
                                         size, rss, share, text, datastack
   else if (masterproc) then
      write(iulog,*)'print_memusage: get_memusage returns -1'
   end if

   return
end subroutine print_memusage
