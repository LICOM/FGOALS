!============================================================================
! offlie radiation driver
!
! Francis Vitt -- Created 15 Dec 2009
!============================================================================

program radiation_driver

  use rad_driver_routines, only: rad_driver_init, rad_driver_final

  use shr_kind_mod,    only: r8 => shr_kind_r8, cl=>shr_kind_cl, cs=>SHR_KIND_CS
  use perf_mod
  use radiation,       only: radiation_tend
  use phys_buffer,     only: pbuf
  use ppgrid,          only: pcols, begchunk, endchunk
  use rad_data_input,  only: ntimes, get_rad_data_input
  use rad_data_input,  only: dates, secs
  use cam_history,     only: wshist, wrapup  
  use physics_types,   only: physics_state, physics_ptend
  use camsrfexch_types,only: cam_in_t, cam_out_t
  use time_manager,    only: timemgr_set_date_time

  use spmd_utils,      only: masterproc, iam
  use cam_logfile,     only: iulog

  use shr_mpi_mod,     only: shr_mpi_barrier
  use shr_mpi_mod,     only: shr_mpi_min, shr_mpi_max
  use shr_mem_mod,     only: shr_mem_init, shr_mem_getusage

  use shr_mpi_mod,     only: shr_mpi_init, shr_mpi_chkerr
  use solar_data,      only: solar_data_advance

  implicit none
#include <mpif.h>

  type(cam_out_t),     pointer :: cam_out(:)       ! Output from CAM to surface
  type(cam_in_t) ,     pointer :: cam_in(:)        ! Merged input state to CAM
  type(physics_state), pointer :: phys_state(:)
  real(r8),            pointer :: landm(:,:)       ! land fraction ramp
  type(physics_ptend)          :: ptend

  real(r8) :: fsns(pcols)      ! Surface solar absorbed flux
  real(r8) :: fsnt(pcols)      ! Net column abs solar flux at model top
  real(r8) :: flns(pcols)      ! Srf longwave cooling (up-down) flux
  real(r8) :: flnt(pcols)      ! Net outgoing lw flux at model top
  real(r8) :: fsds(pcols)      ! Surface solar down flux
  real(r8) :: net_flx(pcols)

  integer :: i,c, irec
  character(len=128) :: buffer
  character(len=128) :: infile
  character(len=128) :: outfile
  logical :: rstwr, nlend

  character(*), parameter :: NLFileName = "drv_in"  ! input namelist filename

  !----------------------------------------------------------------------------
  ! timing vars
  !----------------------------------------------------------------------------
  real(r8)      :: Time_begin        ! Start time
  real(r8)      :: Time_end          ! Ending time

  !----------------------------------------------------------------------------
  ! memory monitoring
  !----------------------------------------------------------------------------
  real(r8) :: msize,msize0,msize1     ! memory size (high water)
  real(r8) :: mrss ,mrss0 ,mrss1      ! resident size (current memory use)

  !----------------------------------------------------------------------------
  ! formats
  !----------------------------------------------------------------------------
  character(*), parameter :: subname = '(rad_driver)'
  character(*), parameter :: FormatR = '(A,": =============== ", A31,F9.3,1x,  " ===============")'

105 format( A, f10.2, A, f10.2, A, i5, A)

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------
   
   if (masterproc) write(iulog,*) 'RADIATION DRIVER BEGIN'

   call shr_mpi_init()
   Time_begin = mpi_wtime()

   call rad_driver_init(cam_in, cam_out, phys_state, landm)

   call t_initf(NLFileName)

   call t_startf('RADIATION_DRIVER')

   !-----------------------------------------------------------------------------
   ! Memory test
   !-----------------------------------------------------------------------------
   call shr_mem_init(prt=.true.)

   !-----------------------------------------------------------------------------
   ! loop through the time records
   !-----------------------------------------------------------------------------
   time_loop: do irec = 1,ntimes

      call t_startf('DRIVER_TIME_LOOP')

      call timemgr_set_date_time(dates(irec),secs(irec))
      call solar_data_advance()

      call t_startf('DRIVER_INPUT_DATA')
      call get_rad_data_input( phys_state, pbuf, cam_in, landm, recno=irec )
      call t_stopf('DRIVER_INPUT_DATA')

      call t_startf('DRIVER_CHUNK_LOOP')
      !$OMP PARALLEL DO PRIVATE (c, fsns, fsnt, flns, flnt, fsds, net_flx)
      chunk_loop: do c=begchunk, endchunk

         call radiation_tend(phys_state(c), ptend, pbuf, cam_out(c), cam_in(c), &
              cam_in(c)%landfrac,landm(:,c),cam_in(c)%icefrac, cam_in(c)%snowhland, &
              fsns, fsnt, flns, flnt, fsds, net_flx)

      enddo chunk_loop
      call t_stopf('DRIVER_CHUNK_LOOP')

      nlend = irec == ntimes
      rstwr = nlend

      if ( masterproc ) then
         call shr_mem_getusage(msize,mrss)
         write(iulog,105) ' memory_write: memory = ',msize,&
              ' MB (highwater)    ',mrss,' MB (usage) (pe=',iam,')'
      endif

      call shr_mpi_barrier(MPI_COMM_WORLD)
      call wshist ()
      call shr_mpi_barrier(MPI_COMM_WORLD)
      call wrapup(rstwr, nlend)
      call shr_mpi_barrier(MPI_COMM_WORLD)

      call t_stopf('DRIVER_TIME_LOOP')

   enddo time_loop

   call t_barrierf ('DRIVER_FINAL_BARRIER', MPI_COMM_WORLD)

   Time_end = mpi_wtime()
   call shr_mem_getusage(msize,mrss)

   call shr_mpi_min(msize,msize0,MPI_COMM_WORLD,'driver msize0',all=.true.)
   call shr_mpi_max(msize,msize1,MPI_COMM_WORLD,'driver msize1',all=.true.)
   call shr_mpi_min(mrss,mrss0,MPI_COMM_WORLD,'driver mrss0',all=.true.)
   call shr_mpi_max(mrss,mrss1,MPI_COMM_WORLD,'driver mrss1',all=.true.)

   if ( masterproc ) then
      write(iulog,FormatR) subname,'compute time (hrs)          = ', (Time_end-Time_begin)/3600._r8
      write(iulog,FormatR) subname,' pes min memory highwater  (MB)  = ',msize0
      write(iulog,FormatR) subname,' pes max memory highwater  (MB)  = ',msize1
      write(iulog,FormatR) subname,' pes min memory last usage (MB)  = ',mrss0
      write(iulog,FormatR) subname,' pes max memory last usage (MB)  = ',mrss1
   endif

   if (masterproc) write(iulog,*) 'RADIATION DRIVER END'

   call t_stopf('RADIATION_DRIVER')
   call t_prf('./rad_drv_timing', MPI_COMM_WORLD )
   call t_finalizef()

   deallocate ( landm, phys_state, cam_in, cam_out )
   call rad_driver_final()
   close(iulog)

end program radiation_driver
