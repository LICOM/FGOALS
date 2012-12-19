!
!BOP
!
! !MODULE: ice_diagnostics - diagnostic information output during run
!
! !DESCRIPTION:
!
! Diagnostic information output during run
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL 
!          Bruce P. Briegleb, NCAR
!
! !INTERFACE:
!
      module ice_diagnostics
!
! !USES:
!
      use ice_constants
!
!EOP
!
      implicit none

      ! for water and heat budgets
      real (kind=dbl_kind) ::
     &   totmn   ! total ice/snow water mass (nh)
     &,  totms   ! total ice/snow water mass (sh)
     &,  toten   ! total ice/snow enthalpy (J)
     &,  totes   ! total ice/snow enthalpy (J)

      ! printing info for routine print_state
      integer, parameter ::
     &  check_step = 100000
     &, ip = 51
     &, jp = 37
     &, mtask = 2
!
      contains
!
!BOP
!
! !IROUTINE: print_state - print ice state for specified grid point
!
! !INTERFACE:
!
      subroutine print_state( plabel,i,c, aice, sicthk, snowhice
     &,          frzmlt, ice_in )
!
! !DESCRIPTION:
!
! This routine is useful for debugging call to it should be inserted 
! in the form (after thermo, for example)
!
!      do j=jlo,jhi    \\
!      -do i=ilo,ihi   \\ 
!      ---call to\_column(i,j)    \\
!      ---plabel = 'post thermo' \\
!      ---if (istep1.ge.check\_step.and.i.eq.ip.and.j.eq.jp  \\
!      ------.and.my\_task.eq.mtask)  \\
!      ---call print\_state(plabel,i,j) \\
!      -enddo  \\
!      enddo   \\
!
! 'use ice\_diagnostics' may need to be inserted also, and \\
! 'use ice\_calendar' if it is not already being used
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_kinds_mod
      use ice_types,        only: ice_in_t
      use phys_grid,        only: get_lat_p, get_lon_p
      use ice_time_manager, only: get_nstep
      use ppgrid,           only: begchunk, endchunk
!
! !INPUT/OUTPUT PARAMETERS:
!
      character (len=20), intent(in) :: plabel
      integer, intent(in) :: i,c
      real(dbl_kind), intent(in) :: aice            ! Ice fraction
      real(dbl_kind), intent(in) :: sicthk          ! Sea-ice thickness
      real(dbl_kind), intent(in) :: snowhice        ! snow height over ice
      real(dbl_kind), intent(in) :: frzmlt          ! Freeze/melt potential
      type(ice_in_t), intent(in) :: ice_in(begchunk:endchunk)
!
!EOP
!
      real (kind=dbl_kind) :: edebug

      integer :: nc, layer
      integer :: nunit
      integer :: lon, lat
      real (kind=dbl_kind) :: z

      lat = get_lat_p (c,i)
      lon = get_lon_p (c,i)

      nunit = 6

      write (nunit,*) plabel,' nstep = ', get_nstep ()
      write (nunit,*) 'lon,lat=',lon,lat
      write (nunit,*) 'aice ',aice
!      write (nunit,*) 'tmask =',tmask(i,j)

      write (nunit,*) 'nc=',1
      write (nunit,*) 'ain ',aice
      write (nunit,*) 'vin ',aice*sicthk
      write (nunit,*) 'hin ',sicthk
      write (nunit,*) 'vsn ',aice*snowhice
      write (nunit,*) 'hsn ',snowhice
!      write (nunit,*) 'Tsfn',Tsfn(nc)
      write (nunit,*) ' '

!      edebug = c0
!      do layer = 1,ni
!        write (nunit,*) 'ein, cat ',nc,' layer ',layer,ein(layer,nc)
!        edebug = edebug+ein(layer,nc)
!      enddo
!      write (nunit,*) 'eice(i,j)',edebug 

      write (nunit,*) ' '
      write (nunit,*) 'atm states and fluxes'
      write (nunit,*) '            zlvl  = ',ice_in%zbot(i)
!     write (nunit,*) '            uatm  = ',uatm(i,c)
!     write (nunit,*) '            vatm  = ',vatm(i,c)
      write (nunit,*) '            potT  = ',ice_in(c)%thbot(i)
      write (nunit,*) '            Tair  = ',ice_in(c)%tbot(i)
      write (nunit,*) '            Qa    = ',ice_in(c)%qbot(i)
!     write (nunit,*) '            rhoa  = ',rhoa(i,c)
      write (nunit,*) '            swvdr = ',ice_in(c)%sols(i)
      write (nunit,*) '            swvdf = ',ice_in(c)%solsd(i)
      write (nunit,*) '            swidr = ',ice_in(c)%soll(i)
      write (nunit,*) '            swidf = ',ice_in(c)%solld(i)
      write (nunit,*) '            Flw   = ',ice_in(c)%flwds(i)
!     write (nunit,*) '            Frain = ',Frain(i,c)
      write (nunit,*) '            Fsnow = ',ice_in(c)%snow(i)
      write (nunit,*) ' '
      write (nunit,*) 'ocn states and fluxes'
      write (nunit,*) '            frzmlt = ',frzmlt
      write (nunit,*) '            sst    = ',ice_in(c)%tsocn(i)
!     write (nunit,*) '            sss    = ',sss(i,c)
      write (nunit,*) '            Tf     = ',-1.8_dbl_kind
      write (nunit,*) ' '
      write (nunit,*) 'srf states and fluxes'
!     write (nunit,*) '            Tref      = ',Tref(i,c)
!     write (nunit,*) '            Fsensible = ',Fsensible(i,c)
!     write (nunit,*) '            Flatent   = ',Flatent(i,c)
!     write (nunit,*) '            evap      = ',evap(i,c)
!     write (nunit,*) '            Flwout    = ',Flwout(i,c)
      write (nunit,*) ' '

      end subroutine print_state
!
      end module ice_diagnostics
!
