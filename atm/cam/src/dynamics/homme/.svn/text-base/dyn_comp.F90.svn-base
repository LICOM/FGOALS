!BOP
!
! !MODULE: dyn_comp --- Dynamical Core Component
!
! !INTERFACE:

Module dyn_comp
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use domain_mod, only : domain1d_t
  use element_mod, only : element_t, elem_state_t
  use time_mod, only : TimeLevel_t
  use hybvcoord_mod, only : hvcoord_t
  use hybrid_mod, only : hybrid_t
  use perf_mod, only: t_startf, t_stopf
  use cam_logfile, only : iulog
  use time_manager,     only: is_first_step
  use spmd_utils,  only : iam, npes
  implicit none
  private


  ! !PUBLIC MEMBER FUNCTIONS:
  public dyn_init1, dyn_init2, dyn_run, dyn_final

  ! !PUBLIC DATA MEMBERS:
  public dyn_import_t, dyn_export_t


  type (TimeLevel_t)   , public :: TimeLevel     ! main time level struct (used by tracers)

!  type (elem_state_t), save, target :: dyn_state

  type dyn_import_t
     type (element_t), pointer :: elem(:)
  end type dyn_import_t

  type dyn_export_t
     type (element_t), pointer :: elem(:)
  end type dyn_export_t
  type (hvcoord_t) :: hvcoord
  integer, parameter  ::  DYN_RUN_SUCCESS           = 0
  integer, parameter  ::  DYN_RUN_FAILURE           = -1

  ! !DESCRIPTION: This module implements the HOMME Dynamical Core as
  !               an ESMF gridded component.  It is specific to HOMME
  !               and does not use ESMF.
  !
  ! \paragraph{Overview}
  !
  !   This module contains an ESMF wrapper for the Homme
  !   Dynamical Core used in the Community Atmospheric Model. 
  !
  ! !REVISION HISTORY:
  !
  !  JPE  06.05.31:  created
  !
  !EOP
  !----------------------------------------------------------------------
  !BOC

  ! Enumeration of DYNAMICS_IN_COUPLINGS


  logical, parameter         :: DEBUG = .true.

  real(r8), parameter        :: ONE    = 1.0D0

  character(*), parameter, public :: MODULE_NAME = "dyn_comp"
  character(*), parameter, public :: VERSION     = "$Id$" 
  type (domain1d_t), pointer, public :: dom_mt(:)      


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-------------------------------------------------------------------------
  !BOP
  ! !ROUTINE:  dyn_init --- Initialize the dynamical core
  !
  ! !INTERFACE:
  subroutine dyn_init1(NLFileName, par, dyn_in, dyn_out)

    ! USES:
    use hycoef,           only: hycoef_init
    use dimensions_mod,   only: globaluniquecols, nelem, nelemd, nelemdmax

    use prim_driver_mod,  only: prim_init1
    use pmgrid,           only: dyndecomp_set
    use dyn_grid,         only : elem, set_horiz_grid_cnt_d, get_dyn_grid_parm, set_horiz_grid_cnt_d
    use thread_mod,       only: nthreads
    use parallel_mod,     only: parallel_t
    use rgrid, only : fullgrid
    use time_mod,         only: TimeLevel_Init
    use spmd_utils,      only: mpi_integer, mpicom, mpi_logical
    use interpolate_mod,  only : interpolate_analysis
    use native_mapping, only : create_native_mapping_files, native_mapping_readnl
    implicit none


    ! PARAMETERS:
    character(len=*), intent(in) :: NLFileName
    type (parallel_t), intent(in) :: par
    type (dyn_import_t), intent(OUT)  :: dyn_in
    type (dyn_export_t), intent(OUT)  :: dyn_out

#ifdef _OPENMP    
    integer omp_get_num_threads
#endif
    integer neltmp(3)
    logical nellogtmp(7)


     !
     ! Initialize dynamics grid
     !

    call hycoef_init
    fullgrid=.true.

#ifdef _OPENMP    
!   Set by driver
!$omp parallel
    nthreads = omp_get_num_threads()
!$omp end parallel
    if(par%masterproc) then
       write(iulog,*) " "
       write(iulog,*) "dyn_init1: nthreads= ",nthreads
       write(iulog,*) " "
    endif
#else
    nthreads = 1
    if(par%masterproc) then
       write(iulog,*) " "
       write(iulog,*) "dyn_init1: openmp not activated"
       write(iulog,*) " "
    endif
#endif
    if(iam < par%nprocs) then
       call prim_init1(elem,par,dom_mt,TimeLevel)

       dyn_in%elem => elem
       dyn_out%elem => elem
    
       call set_horiz_grid_cnt_d(GlobalUniqueCols)


       neltmp(1) = nelemdmax
       neltmp(2) = nelem
       neltmp(3) = get_dyn_grid_parm('plon')
       nellogtmp(1:7) = interpolate_analysis(1:7)
    else
       nelemd = 0
       neltmp(1) = 0
       neltmp(2) = 0
       neltmp(3) = 0
       nellogtmp(1:7) = .true.
    endif

    dyndecomp_set = .true.



    if (par%nprocs .lt. npes) then
! Broadcast quantities to auxiliary processes
       call mpibcast(neltmp, 3, mpi_integer, 0, mpicom)
       call mpibcast(nellogtmp, 7, mpi_logical, 0, mpicom)
       if (iam .ge. par%nprocs) then
          nelemdmax = neltmp(1)
          nelem     = neltmp(2)
          call set_horiz_grid_cnt_d(neltmp(3))
          interpolate_analysis(1:7) = nellogtmp(1:7)
       endif
    endif


    !
    ! This subroutine creates mapping files using homme basis functions if requested
    !
    call native_mapping_readnl(NLFileName)
    call create_native_mapping_files( par, elem,'native')
    call create_native_mapping_files( par, elem,'bilin')

  end subroutine dyn_init1


  subroutine dyn_init2(par, elem)
    use dimensions_mod,   only: nlev, nelemd
    use prim_driver_mod,  only: prim_init2, prim_run
    use prim_si_ref_mod,  only: prim_set_mass
    use hybrid_mod,       only: hybrid_create
    use hycoef,           only: hyam, hybm, hyai, hybi, ps0
    use parallel_mod,     only: parallel_t
    use time_manager,     only: dtime,get_nstep   ! physics timestep
    use time_mod,         only: homme_nsplit, tstep, time_at
    use control_mod,      only: moisture, runtype, qsplit
    use thread_mod,       only: nthreads
    use thread_mod,       only: omp_get_thread_num
    use cam_control_mod,  only: aqua_planet, ideal_phys, adiabatic
    use comsrf,           only: landm, sgh, sgh30
    use nctopo_util_mod,  only: nctopo_util_driver


    type(element_t), intent(inout) :: elem(:)
    type(parallel_t), intent(in) :: par
    integer :: ithr, nets, nete, ie, k
    real(r8), parameter :: Tinit=300.0_r8
    real(r8) :: dyn_ps0
    type(hybrid_t) :: hybrid

    !
    !  Note: dtime = progress made in one timestep.  value in namelist
    !        dtime = the frequency at which physics is called
    !        tstep = the dynamics timestep:  
    !
    !    Leapfrog looks like:   u(3) = u(1) + 2*tstep*u(2)   
    !    u(1) = time-tstep 
    !    u(2) = time
    !    u(3) = time+tstep 
    !  
    !    Physics looks like:    u(3) = u(1) + dt_phys*PHYSICS(U(1))
    !  
    ! so with homme_nsplit=1:    dtime=tstep   dt_phys=2*tstep
    !
    ! In general:  dtime=homme_nsplit*tstep,  dt_phys=homme_nsplit*tstep + tstep
    ! 

    dyn_ps0=ps0/100.D0
    hvcoord%hyam=hyam
    hvcoord%hyai=hyai
    hvcoord%hybm=hybm
    hvcoord%hybi=hybi
    hvcoord%ps0=dyn_ps0  
    do k=1,nlev
       hvcoord%hybd(k) = hvcoord%hybi(k+1) - hvcoord%hybi(k)
    end do
    if(iam < par%nprocs) then

       !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ie,ithr,nets,nete,hybrid)
       ithr=omp_get_thread_num()
       nets=dom_mt(ithr)%start
       nete=dom_mt(ithr)%end
       hybrid = hybrid_create(par,ithr,NThreads)

       tstep = dtime/float(homme_nsplit*qsplit)

       moisture='moist'

       if(adiabatic) then
          moisture='dry'
          if(runtype == 0) then
             do ie=nets,nete
                elem(ie)%state%q(:,:,:,:,:)=0.0D0
                elem(ie)%derived%fq(:,:,:,:,:)=0.0D0
             end do
          end if
       else if(ideal_phys) then
          moisture='dry'
          if(runtype == 0) then
             do ie=nets,nete
                elem(ie)%state%lnps(:,:,:) =LOG(dyn_ps0)

                elem(ie)%state%ps_v(:,:,:) =dyn_ps0

                elem(ie)%state%phis(:,:)=0.0D0

                elem(ie)%state%T(:,:,:,:) =Tinit

                elem(ie)%state%v(:,:,:,:,:) =0.0D0

                elem(ie)%state%q(:,:,:,:,:)=0.0D0

             end do
          end if
       else if(aqua_planet .and. runtype==0)  then
          do ie=nets,nete
             !          elem(ie)%state%lnps(:,:,:) =LOG(dyn_ps0)
             !          elem(ie)%state%ps_v(:,:,:) =dyn_ps0
             elem(ie)%state%phis(:,:)=0.0D0
          end do
          if(allocated(landm)) landm=0.0_r8
          if(allocated(sgh)) sgh=0.0_r8
          if(allocated(sgh30)) sgh30=0.0_r8
       end if

       do ie=nets,nete
          elem(ie)%derived%FM=0.0D0
          elem(ie)%derived%FT=0.0D0
          elem(ie)%derived%FQ=0.0D0
       end do

       ! initial homme (subcycled) nstep
       TimeLevel%nstep = get_nstep()*homme_nsplit*qsplit

       ! scale PS to achieve prescribed dry mass
       if (runtype == 1) then
          ! exact restart
          TimeLevel%nstep0=TimeLevel%nstep+1
       else
          ! new run, scale mass to value given in namelist, if needed
          call prim_set_mass(elem, TimeLevel,hybrid,hvcoord,nets,nete)  
          TimeLevel%nstep0=2   ! This will be the first full leapfrog step
       endif
       call prim_init2(elem,hybrid,nets,nete, TimeLevel, hvcoord)
       !
       ! This subroutine is used to create nc_topo files, if requested
       ! 
       call nctopo_util_driver(elem,hybrid,nets,nete)
       !$OMP END PARALLEL 
    end if


    
  end subroutine dyn_init2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------
  !BOP
  ! !ROUTINE:  RUN --- Driver for the 
  !
  ! !INTERFACE:
  subroutine dyn_run( dyn_state, rc )

    ! !USES:
    use parallel_mod,     only : par
    use prim_driver_mod,  only: prim_run, prim_run_subcycle
    use dimensions_mod,   only : nlev
    use thread_mod,       only: omp_get_thread_num, nthreads
    use time_mod,         only: homme_nsplit, tstep
    use control_mod,      only: tstep_type
    use hybrid_mod,       only: hybrid_create
!    use perf_mod, only : t_startf, t_stopf
    implicit none


    type (dyn_export_t), intent(inout)       :: dyn_state   !  container
    type(hybrid_t) :: hybrid

    integer, intent(out)               :: rc      ! Return code
    integer ::  n
    integer :: nets, nete, ithr
    integer :: ie
    real(r8) :: tstep_tmp

    ! !DESCRIPTION:
    !
    if(iam < par%nprocs) then
       !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid,n)
       ithr=omp_get_thread_num()
       nets=dom_mt(ithr)%start
       nete=dom_mt(ithr)%end
       hybrid = hybrid_create(par,ithr,NThreads)

       do n=1,homme_nsplit
          if (tstep_type==1) then
             ! forward-in-time RK, with subcycling
             call prim_run_subcycle(dyn_state%elem, hybrid,nets,nete, tstep, TimeLevel, hvcoord)
          else
             ! leapfrog
             call prim_run(dyn_state%elem, hybrid,nets,nete, tstep, TimeLevel, hvcoord, "leapfrog")
          endif
       end do


       !$OMP END PARALLEL
    end if
    rc = DYN_RUN_SUCCESS

    !EOC
  end subroutine dyn_run
  !-----------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dyn_final(DYN_STATE, RESTART_FILE)

    type (elem_state_t), target     :: DYN_STATE
    character(LEN=*)   , intent(IN) :: RESTART_FILE



  end subroutine dyn_final

end module dyn_comp



