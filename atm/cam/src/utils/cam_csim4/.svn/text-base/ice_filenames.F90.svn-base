module ice_filenames

!----------------------------------------------------------------------- 
!
! DESCRIPTION
! Module and methods to handle filenames needed for the model. This 
! includes input filenames, and most output filenames that the model
! uses. All filenames that the model uses will use methods or data
! constructed by this module. In some cases (such as the cam_history module)
! other modules or routines will store the actual filenames used, but
! this module is used to determine the names.
!
!----------------------------------------------------------------------- 

   use ice_time_manager, only: get_curr_date, get_prev_date
   use shr_kind_mod,     only: shr_kind_cs, shr_kind_cl
   use shr_sys_mod,      only: shr_sys_abort
   use cam_logfile,      only: iulog

   use filenames,        only: caseid    !cam specific use

   implicit none

   integer, parameter :: nlen = shr_kind_cl                ! String length

CONTAINS

  character(len=nlen) function interpret_filename_spec( filename_spec, number, prev, case, &
       yr_spec, mon_spec, day_spec, sec_spec )

  !----------------------------------------------------------------------- 
  !
  ! !DESCRIPTION: Create a filename from a filename specifyer. The 
  ! filename specifyer includes codes for setting things such as the
  ! year, month, day, seconds in day, caseid, and tape number. This
  ! routine is private to filenames.F90
  !
  ! Interpret filename specifyer string with: 
  !
  !      %c for case, 
  !      %t for optional number argument sent into function
  !      %y for year
  !      %m for month
  !      %d for day
  !      %s for second
  !      %% for the "%" character
  !
  ! If the filename specifyer has spaces " ", they will be trimmed out
  ! of the resulting filename.
  !
  !----------------------------------------------------------------------- 
  !
  ! Arguments
  !
  character(len=*), intent(in)           :: filename_spec  ! Filename specifier to use
  integer         , intent(in), optional :: number         ! Number to use for %t field
  logical         , intent(in), optional :: prev           ! If should label with previous time-step
  character(len=*), intent(in), optional :: case           ! Optional casename
  integer         , intent(in), optional :: yr_spec        ! Simulation year
  integer         , intent(in), optional :: mon_spec       ! Simulation month
  integer         , intent(in), optional :: day_spec       ! Simulation day
  integer         , intent(in), optional :: sec_spec       ! Seconds into current simulation day
  !
  ! Local variables
  !
  integer :: year  ! Simulation year
  integer :: month ! Simulation month
  integer :: day   ! Simulation day
  integer :: ncsec ! Seconds into current simulation day
  character(len=nlen) :: string    ! Temporary character string 
  character(len=nlen) :: format    ! Format character string 
  integer :: i, n                  ! Loop variables
  logical :: previous              ! If should label with previous time-step
  logical :: done
  !----------------------------------------------------------------------- 

  if ( len_trim(filename_spec) == 0 )then
     call shr_sys_abort ('INTERPRET_FILENAME_SPEC: filename specifier is empty')
  end if
  if ( index(trim(filename_spec)," ") /= 0 )then
     call shr_sys_abort ('INTERPRET_FILENAME_SPEC: filename specifier can not contain a space:'//trim(filename_spec))
  end if
  !
  ! Determine year, month, day and sec to put in filename
  !
  if (present(yr_spec) .and. present(mon_spec) .and. present(day_spec) .and. present(sec_spec)) then
     year  = yr_spec
     month = mon_spec
     day   = day_spec
     ncsec = sec_spec
  else
     if ( .not. present(prev) ) then
        previous = .false.
     else
        previous = prev
     end if
     if ( previous ) then
        call get_prev_date(year, month, day, ncsec)
     else
        call get_curr_date(year, month, day, ncsec)
     end if
  end if
  !
  ! Go through each character in the filename specifyer and interpret if special string
  !
  i = 1
  interpret_filename_spec = ''
  do while ( i <= len_trim(filename_spec) )
     !
     ! If following is an expansion string
     !
     if ( filename_spec(i:i) == "%" )then
        i = i + 1
        select case( filename_spec(i:i) )
        case( 'c' )   ! caseid
           if ( present(case) )then
              string = trim(case)
           else
              string = trim(caseid)
           end if
        case( 't' )   ! number
           if ( .not. present(number) )then
              write(iulog,*) 'INTERPRET_FILENAME_SPEC: number needed in filename_spec' &
                   , ', but not provided to subroutine'
              write(iulog,*) 'filename_spec = ', filename_spec
              call shr_sys_abort
           end if
           if (      number > 999 ) then
              format = '(i4.4)'
              if ( number > 9999 ) then
                 write(iulog,*) 'INTERPRET_FILENAME_SPEC: number is too large: ', number
                 call shr_sys_abort
              end if
           else if ( number > 99  ) then
              format = '(i3.3)'
           else if ( number > 9   ) then
              format = '(i2.2)'
           else
              format = '(i1.1)'
           end if
           write(string,format) number
        case( 'y' )   ! year
           if ( year > 99999   ) then
              format = '(i6.6)'
           else if ( year > 9999    ) then
              format = '(i5.5)'
           else
              format = '(i4.4)'
           end if
           write(string,format) year
        case( 'm' )   ! month
           write(string,'(i2.2)') month
        case( 'd' )   ! day
           write(string,'(i2.2)') day
        case( 's' )   ! second
           write(string,'(i5.5)') ncsec
        case( '%' )   ! percent character
           string = "%"
        case default
           call shr_sys_abort ('INTERPRET_FILENAME_SPEC: Invalid expansion character: '//filename_spec(i:i))
        end select
        !
        ! Otherwise take normal text up to the next "%" character
        !
     else
        n = index( filename_spec(i:), "%" )
        if ( n == 0 ) n = len_trim( filename_spec(i:) ) + 1
        if ( n == 0 ) exit 
        string = filename_spec(i:n+i-2)
        i = n + i - 2
     end if
     if ( len_trim(interpret_filename_spec) == 0 )then
        interpret_filename_spec = trim(string)
     else
        if ( (len_trim(interpret_filename_spec)+len_trim(string)) >= nlen )then
           call shr_sys_abort ('INTERPRET_FILENAME_SPEC: Resultant filename too long')
        end if
        interpret_filename_spec = trim(interpret_filename_spec) // trim(string)
     end if
     i = i + 1

  end do
  if ( len_trim(interpret_filename_spec) == 0 )then
     call shr_sys_abort('INTERPRET_FILENAME_SPEC: Resulting filename is empty')
  end if

end function interpret_filename_spec

end module ice_filenames
