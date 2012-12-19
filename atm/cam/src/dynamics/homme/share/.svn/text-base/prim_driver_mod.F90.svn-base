#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"file: ",__FILE__," line: ",__LINE__," ithr: ",hybrid%ithr
#define _DBG_ 
module prim_driver_mod
  use kinds, only : real_kind, iulog
  use cg_mod, only : cg_t
  use hybrid_mod, only : hybrid_t
  use quadrature_mod, only : quadrature_t
#ifndef CAM
  use column_model_mod, only : ColumnModel_t
  use prim_restart_mod, only : initrestartfile
  use restart_io_mod , only : RestFile,readrestart
#endif
  use prim_si_ref_mod, only : ref_state_t
  use solver_mod, only : blkjac_t
  use filter_mod, only : filter_t
  use derivative_mod, only : derivative_t
  use reduction_mod, only : reductionbuffer_ordered_1d_t, red_min, red_max, &
         red_sum, red_flops, initreductionbuffer,parallelsum
  implicit none
  private
  public :: prim_init1, prim_init2 , prim_run, prim_run_subcycle, prim_finalize, leapfrog_bootstrap
  public :: smooth_topo_datasets

  type (cg_t), allocatable  :: cg(:)              ! conjugate gradient struct (nthreads)

  type (quadrature_t)   :: gv,gp           ! quadratures on velocity and pressure grids
#ifndef CAM
  type (ColumnModel_t), allocatable :: cm(:) ! (nthreads)
#endif
  type (ref_state_t)    :: refstate        ! semi-implicit vertical reference state
  type (blkjac_t),allocatable  :: blkjac(:)  ! (nets:nete)
  type (filter_t)       :: flt             ! Filter struct for v and p grid
  type (filter_t)       :: flt_advection   ! Filter struct for v grid for advection only
  type (derivative_t), allocatable   :: deriv(:) ! derivative struct (nthreads)
  real*8  :: tot_iter
  type (ReductionBuffer_ordered_1d_t) :: red   ! reduction buffer               (shared)

contains
  subroutine prim_init1(elem, par, dom_mt, Tl)
    ! --------------------------------
    use element_mod, only : element_t
    use thread_mod, only : nthreads, omp_get_thread_num, omp_set_num_threads
    ! --------------------------------
    use control_mod, only : runtype, restartfreq, filter_counter, integration, topology, &
         partmethod, while_iter
    ! --------------------------------
    use prim_state_mod, only : prim_printstate_init
    ! --------------------------------
    use namelist_mod, only : readnl
    ! --------------------------------
    use dimensions_mod, only : nv, np, nlev, nelem, nelemd, ne, nelemdmax, GlobalUniqueCols
    ! -------------------------------- 
    use time_mod, only : nmax, time_at, timelevel_init, timelevel_t
    ! --------------------------------
    use quadrature_mod, only : test_gauss, test_gausslobatto
    ! --------------------------------
    use element_mod, only : element_t
    ! --------------------------------
    use mass_matrix_mod, only : mass_matrix
    ! --------------------------------
    use cube_mod,  only : cubeedgecount , cubeelemcount, cubetopology, &
         cube_init_atomic, rotation_init_atomic
    ! --------------------------------
    use metagraph_mod, only : metavertex_t, metaedge_t, localelemcount, initmetagraph
    ! --------------------------------
    use gridgraph_mod, only : gridvertex_t, gridedge_t
    ! --------------------------------
    use schedule_mod, only : schedule, genEdgeSched
    ! --------------------------------
    use prim_advection_mod, only: prim_advec_init
    ! --------------------------------    
    use prim_advance_mod, only: prim_advance_init
    ! --------------------------------    
    use diffusion_mod, only      : diffusion_init
    ! --------------------------------    
    use parallel_mod, only : iam, parallel_t, syncmp, abortmp, global_shared_buf
#ifdef _MPI
    use parallel_mod, only : mpiinteger_t, mpireal_t, mpi_max, mpi_sum, haltmp
#endif
    ! --------------------------------
    use metis_mod, only : genmetispart
    ! --------------------------------
#ifdef TESTGRID
    use checksum_mod, only : testchecksum
#endif
    ! --------------------------------
    use spacecurve_mod, only : genspacepart
    ! --------------------------------
    use dof_mod, only : global_dof, CreateUniqueIndex, SetElemOffset
    ! --------------------------------
    use params_mod, only : SFCURVE
    ! --------------------------------
!    use surfaces_mod, only: InitControlVolumes1, InitControlVolumes2, &
!         VerifVolumes
    ! --------------------------------
    use hybrid_mod, only : hybrid_t, hybrid_create
    ! --------------------------------
    use domain_mod, only : domain1d_t, decompose
    ! --------------------------------
    use physical_constants, only : dd_pi
    ! --------------------------------
#ifdef CAM
    use repro_sum_mod, only: repro_sum
#endif
    implicit none
    type (element_t), pointer :: elem(:)

    type (parallel_t), intent(in) :: par
    type (domain1d_t), pointer :: dom_mt(:)      
    type (timelevel_t), intent(out) :: Tl
    ! Local Variables
    type (hybrid_t) :: hybrid

    type (GridVertex_t), target,allocatable :: GridVertex(:)
    type (GridEdge_t),   target,allocatable :: Gridedge(:)
    type (MetaVertex_t), target,allocatable :: MetaVertex(:)
    type (MetaEdge_t),   target,allocatable :: MetaEdge(:)

    integer :: ii,ie, ith
    integer :: nets, nete
    integer :: nelem_edge,nedge
    integer :: nstep
    integer :: nlyr
    integer :: iMv
    integer :: err, ierr

    real(kind=real_kind), allocatable :: aratio(:,:)
    real(kind=real_kind) :: area(1)
    character(len=80) rot_type   ! cube edge rotation type

    integer  :: i
    integer,allocatable :: TailPartition(:)
    integer,allocatable :: HeadPartition(:)

    integer total_nelem
    real approx_elements_per_task

    ! =====================================
    ! Read in model control information
    ! =====================================
    ! cam readnl is called in spmd_dyn (needed prior to mpi_init)
#ifndef CAM
    call readnl(par)
    total_nelem      = CubeElemCount(ne)
    approx_elements_per_task = total_nelem/par%nprocs
    if  (approx_elements_per_task < 1.0D0) then
       if(par%masterproc) print *,"number of elements=", total_nelem
       if(par%masterproc) print *,"number of procs=", par%nprocs
       call abortmp('There is not enough paralellism in the job, that is, there is less than one elements per task.')
    end if
#endif
    ! ====================================
    ! Set cube edge rotation type for model
    ! unnecessary complication here: all should
    ! be on the same footing. RDL
    ! =====================================
    rot_type="contravariant"

    if (par%masterproc) then

       ! =============================================
       ! Compute total simulated time...
       ! =============================================

       write(iulog,*)""
       write(iulog,*)" total simulated time = ",Time_at(nmax)
       write(iulog,*)""

       ! =============================================
       ! Perform Gauss/Gauss Lobatto tests...
       ! =============================================

       call test_gauss(np)
       call test_gausslobatto(nv)

    end if

    ! ===============================================================
    ! Allocate and initialize the graph (array of GridVertex_t types)
    ! ===============================================================

    if (topology=="cube") then

       if (par%masterproc) then
          write(iulog,*)"creating cube topology..."
       end if

       nelem      = CubeElemCount(ne)
       nelem_edge = CubeEdgeCount(ne) 
       allocate(GridVertex(nelem))
       allocate(GridEdge(nelem_edge))

       call CubeTopology(GridEdge,GridVertex)

       if(par%masterproc)       write(iulog,*)"...done."
    end if


    !debug  call PrintGridVertex(GridVertex)

    if(par%masterproc) write(iulog,*)"partitioning graph..."

    if(partmethod .eq. SFCURVE) then 
       call genspacepart(GridEdge,GridVertex)
    else
       call genmetispart(GridEdge,GridVertex)
    endif

    ! ===========================================================
    ! given partition, count number of local element descriptors
    ! ===========================================================
#ifdef _PREDICT
    allocate(MetaVertex(npart))
    allocate(Schedule(npart))
#ifdef _MPI
    call abortmp("init: PREDICT code branch not supported under MPI")
#else
    allocate(elem(nelem))
#endif

    do iMv = 1,npart
       ! ====================================================
       !  Generate the communication graph
       ! ====================================================
       call initMetaGraph(iMv,MetaVertex(iMv),GridVertex,GridEdge)

       nelemd = LocalElemCount(MetaVertex(iMv))

       ! ====================================================
       !  Generate the communication schedule
       ! ====================================================
       call genEdgeSched(iMv,Schedule(iMv),MetaVertex(iMv))

    enddo

    call PrimMessageStats(nlyr)
    stop
#else
    allocate(MetaVertex(1))
    allocate(Schedule(1))

    nelem_edge=SIZE(GridEdge)

    allocate(TailPartition(nelem_edge))
    allocate(HeadPartition(nelem_edge))
    do i=1,nelem_edge
       TailPartition(i)=GridEdge(i)%tail%partition
       HeadPartition(i)=GridEdge(i)%head%partition
    enddo

    ! ====================================================
    !  Generate the communication graph
    ! ====================================================
    call initMetaGraph(iam,MetaVertex(1),GridVertex,GridEdge)

    nelemd = LocalElemCount(MetaVertex(1))

    if(nelemd .le. 0) then 
       call abortmp('Not yet ready to handle nelemd = 0 yet' )
       stop
    endif

#ifdef _MPI 
    call mpi_allreduce(nelemd,nelemdmax,1, &
         MPIinteger_t,MPI_MAX,par%comm,ierr)
#else
    nelemdmax=nelemd
#endif

    allocate(elem(nelemd),stat=err)

    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================

    call genEdgeSched(elem,iam,Schedule(1),MetaVertex(1))

#endif
    allocate(global_shared_buf(nelemd,5), stat=err)
    !  nlyr=edge3p1%nlyr
    !  call MessageStats(nlyr)
    !  call testchecksum(par,GridEdge)

    ! ========================================================
    ! load graph information into local element descriptors
    ! ========================================================

    !  do ii=1,nelemd
    !     elem(ii)%vertex = MetaVertex(iam)%members(ii)
    !  enddo

    call syncmp(par)

    ! =================================================================
    ! Initialize shared boundary_exchange and reduction buffers
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'init shared boundary_exchange buffers'
    call InitReductionBuffer(red,3*nlev,nthreads)
    call InitReductionBuffer(red_sum,5)
    call InitReductionBuffer(red_max,1)
    call InitReductionBuffer(red_min,1)
    call initReductionBuffer(red_flops,1)


    if (topology=="cube") then
       if(par%masterproc) write(iulog,*) "initializing cube elements..."
       do ie=1,nelemd
          call cube_init_atomic(elem(ie))
       enddo
    end if

    ! =================================================================
    ! Initialize mass_matrix
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'running mass_matrix'
    call mass_matrix(par,elem)
    allocate(aratio(nelemd,1))

    if (topology=="cube") then
       area = 0
       do ie=1,nelemd
          aratio(ie,1) = sum(elem(ie)%mv(:,:)*elem(ie)%metdet(:,:))
       enddo
       ! note: cant use parallelsum() since hybrid NOT YET INITIALIZED
#ifdef CAM
       call repro_sum(aratio, area, nelemd, nelemd, 1, commid=par%comm)
#else
       if(par%masterproc) print *, 'Warning - this is not an task count independent result.'
       do ie=2,nelemd
          aratio(1,1)=aratio(1,1)+aratio(ie,1)
       end do
#ifdef _MPI 
       call MPI_Allreduce(aratio(1,1),area(1),1,MPIreal_t,MPI_SUM,par%comm,ierr)       
#else
       area = aratio(1,1)
#endif
#endif
       area(1) = 4*dd_pi/area(1)  ! ratio correction
       deallocate(aratio)
       if (par%masterproc) &
            write(iulog,'(a,f20.18)') " re-initializing cube elements: area correction=",area(1)
       do ie=1,nelemd
          call cube_init_atomic(elem(ie),area(1))
          call rotation_init_atomic(elem(ie),rot_type)
       enddo
    end if
    if(par%masterproc) write(iulog,*) 're-running mass_matrix'
    call mass_matrix(par,elem)
    ! =================================================================
    ! Run the checksum to verify communication schedule
    ! =================================================================
#ifdef TESTGRID 
    if(par%masterproc)     write(iulog,*) 'running testchecksum ',iam,nelem,nelemd
    call testchecksum(par,GridEdge)
#endif

    ! =================================================================
    ! Determine the global degree of freedome for each gridpoint
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'running global_dof'
    call global_dof(par,elem)

    ! =================================================================
    ! Create Unique Indices
    ! =================================================================

    do ie=1,nelemd
       call CreateUniqueIndex(elem(ie)%GlobalId,elem(ie)%gdofP,elem(ie)%idxP)
    enddo

    call SetElemOffset(par,elem, GlobalUniqueCols)

    do ie=1,nelemd
       elem(ie)%idxV=>elem(ie)%idxP
    end do

    !JMD call PrintDofP(elem)
    !JMD call PrintDofV(elem)




    call prim_printstate_init(par)
    ! Initialize output fields for plotting...

    while_iter = 0
    filter_counter = 0

    ! initialize flux terms to 0

    do ie=1,nelemd
       elem(ie)%derived%FM=0.0
       elem(ie)%derived%FQ=0.0
       elem(ie)%derived%FQps=0.0
       elem(ie)%derived%FT=0.0

       elem(ie)%accum%Qvar=0
       elem(ie)%accum%Qmass=0
       elem(ie)%accum%Q1mass=0
       elem(ie)%accum%mass_added=0
    enddo


    ! ==========================================================
    !  This routines initalizes a Restart file.  This involves:
    !      I)  Setting up the MPI datastructures
    ! ==========================================================
#ifndef CAM
    if(restartfreq > 0 .or. runtype>=1)  then 
       call initRestartFile(elem(1)%state,par,RestFile)
    endif
#endif
    !DBG  write(iulog,*) 'prim_init: after call to initRestartFile'

#ifndef TESTGRID
    deallocate(GridEdge)
    deallocate(GridVertex)
#else
    ! here we need to call a function in gridgraph_mod.F to deallocate
    ! all of GridVertex's EdgeIndex pointers
#endif


    ! =====================================
    ! Begin threaded region...
    ! =====================================
    ! =====================================
    ! Set number of threads...
    ! =====================================
    if(par%masterproc) write(iulog,*) "Main:NThreads=",NThreads
    call omp_set_num_threads(NThreads)

    allocate(dom_mt(0:NThreads-1))
    do ith=0,NThreads-1
       dom_mt(ith)=decompose(1,nelemd,NThreads,ith)
    end do
    ith=0
    nets=1
    nete=nelemd
!    call InitControlVolumes1(nelemd)
#ifndef CAM
    allocate(cm(0:Nthreads-1))
#endif
    allocate(deriv(0:Nthreads-1))
    allocate(cg(0:Nthreads-1))
    call prim_advance_init(integration)
    call Prim_Advec_Init()
    call diffusion_init()



    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ith,nets,nete, hybrid)
    ith=omp_get_thread_num()
    nets=dom_mt(ith)%start
    nete=dom_mt(ith)%end

    hybrid = hybrid_create(par,ith,NThreads)

!    call InitControlVolumes2(elem,hybrid,nets,nete)
!    call VerifVolumes(elem,hybrid,nets,nete)

    !$OMP END PARALLEL 
    ! ===========================================================
    ! initialize the time level of the model...
    ! ===========================================================

    call TimeLevel_init(tl)

    if(par%masterproc) write(iulog,*) 'end of prim_init'
  end subroutine prim_init1

  subroutine prim_init2(elem,hybrid, nets, nete, tl, hvcoord)
    use element_mod, only : element_t
    use parallel_mod, only : parallel_t, haltmp, syncmp, abortmp
    use time_mod, only : timelevel_t, tstep, phys_tscale, timelevel_init, time_at, timelevel_update, nendstep, smooth, homme_nsplit
    use prim_state_mod, only : prim_printstate, prim_diag_scalars
    use filter_mod, only : filter_t, fm_filter_create, taylor_filter_create, &
         fm_transfer, bv_transfer
    use control_mod, only : runtype, integration, filter_mu, filter_mu_advection, test_case, &
         debug_level, vfile_int, filter_freq, filter_freq_advection, &
         transfer_type, vform, vfile_mid, filter_type, kcut_fm, wght_fm, p_bv, &
         s_bv, topology,columnpackage, moisture, precon_method, qsplit, rk_stage_user,&
         TRACERADV_TOTAL_DIVERGENCE, TRACERADV_UGRADQ, tracer_advection_formulation, sub_case
    use quadrature_mod, only : gauss, gausslobatto
    use dimensions_mod, only :  nv, np, nlev, ne, qsize
    use prim_si_ref_mod, only: prim_si_refstate_init, prim_set_mass
    use thread_mod, only : nthreads
    use derivative_mod, only : derivinit
    use global_norms_mod, only : test_global_integral
    use hybvcoord_mod, only : hvcoord_t
#ifdef CAM
#else
    use column_model_mod, only : InitColumnModel
    use held_suarez_mod, only : hs0_init_state
    use baroclinic_inst_mod, only : binst_init_state, jw_baroclinic
    use asp_tests_mod, only : asp_tracer, asp_baroclinic, asp_rossby, asp_mountain, asp_gravity_wave ! _EXTERNAL
    use aquaplanet, only : aquaplanet_init_state
#endif

    type (element_t), intent(inout) :: elem(:)
    type (hybrid_t), intent(in) :: hybrid

    type (TimeLevel_t), intent(inout)    :: tl              ! time level struct
    type (hvcoord_t), intent(inout)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)

    ! ==================================
    ! Local variables
    ! ==================================

    real (kind=real_kind) :: dt              ! "timestep dependent" timestep
    real (kind=real_kind) :: dp        


    real (kind=real_kind) :: ps(nv,nv)       ! surface pressure

    character(len=80)     :: fname
    character(len=8)      :: njusn
    character(len=4)      :: charnum

    real (kind=real_kind) :: Tv(nv)     ! transfer function (velocity grid)
    real (kind=real_kind) :: Tp(np)     ! transfer function (pressure grid)

    integer :: simday
    integer :: i,j,k,ie,iptr,t,q
    integer :: ierr
    integer :: nfrc

    ! ==========================
    ! begin executable code
    ! ==========================
    if (topology == "cube") then
       call test_global_integral(elem, hybrid,nets,nete)  
    end if

    ! ========================================
    ! Initialize velocity and pressure grid
    ! quadrature points...
    ! ========================================

    gv=gausslobatto(nv)

    if (nv==np) then
       gp =gausslobatto(np)
    else
       gp =gauss(np)
    end if

    ! ==================================
    ! Initialize derivative structure
    ! ==================================

    call derivinit(deriv(hybrid%ithr))


    ! ====================================
    ! In the semi-implicit case:
    ! initialize vertical structure and 
    ! related matrices..
    ! ====================================
!$OMP MASTER
    if (integration == "semi_imp") then
       refstate = prim_si_refstate_init(.false.,hybrid%masterthread,hvcoord)
       if (precon_method == "block_jacobi") then
          allocate(blkjac(nets:nete))
       endif
    endif
!$OMP END MASTER
    ! ==========================================
    ! Initialize pressure and velocity grid 
    ! filter matrix...
    ! ==========================================


    if (transfer_type == "bv") then
       Tv    = bv_transfer(p_bv,s_bv,nv)
       Tp    = bv_transfer(p_bv,s_bv,np)
    else if (transfer_type == "fm") then
       Tv    = fm_transfer(kcut_fm,wght_fm,nv)
       Tp    = fm_transfer(kcut_fm,wght_fm,np)
    end if
    if (filter_type == "taylor") then
       flt           = taylor_filter_create(Tp,Tv, filter_mu,gp, gv)
       flt_advection = taylor_filter_create(Tp,Tv, filter_mu_advection,gp, gv)
    else if (filter_type == "fischer") then
       flt           = fm_filter_create(Tp,Tv, filter_mu, gp, gv)
       flt_advection = fm_filter_create(Tp,Tv, filter_mu_advection, gp, gv)
    end if

    if (hybrid%masterthread) then
       if (filter_freq>0 .or. filter_freq_advection>0) then
          write(iulog,*) "transfer function type in preq=",transfer_type
          write(iulog,*) "filter type            in preq=",filter_type
          write(*,'(a,99f10.6)') "dynamics: I-mu + mu*Tv(:) = ",&
               (1-filter_mu)+filter_mu*Tv(:)
          write(*,'(a,99f10.6)') "advection: I-mu + mu*Tv(:) = ",&
               (1-filter_mu_advection)+filter_mu_advection*Tv(:)
       endif
    endif

    !$OMP BARRIER
    if (hybrid%ithr==0) then
       call syncmp(hybrid%par)
    end if
    !$OMP BARRIER

    if (topology /= "cube") then
       call abortmp('Error: only cube topology supported for primaitve equations') 
    endif




#ifndef CAM
    ! =================================
    ! HOMME stand alone initialization
    ! =================================
    tl%nstep0=2   ! This will be the first full leapfrog step
    call InitColumnModel(elem, cm(hybrid%ithr),hvcoord,tl,nets,nete,runtype)

    if(runtype >= 1) then 
       ! ===========================================================
       ! runtype==1   Exact Restart 
       ! runtype==2   Initial run, but take inital condition from Restart file
       ! ===========================================================
       if (hybrid%masterthread) then
          write(iulog,*) 'runtype: RESTART of primitive equations'
       end if
       
       if (test_case(1:10) == "aquaplanet") then
          if (hybrid%masterthread) then
             write(iulog,*)  'Initializing aqua planet with MJO-type forcing'
          end if
          if(moisture.eq."dry") then
             call binst_init_state(elem, hybrid,nets,nete,hvcoord)
          end if
          call aquaplanet_init_state(elem, hybrid,hvcoord,nets,nete,integration)
       end if
       
       call ReadRestart(elem,hybrid%ithr,nets,nete,tl)

       ! scale PS to achieve prescribed dry mass
       if (runtype /= 1) &
            call prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)  
       
       tl%nstep0=tl%nstep+1            ! restart run: first step = first first full leapfrog step
       
       
       if (runtype==2) then
          ! branch run
          ! reset time counters to zero since timestep may have changed
          nEndStep = nEndStep-tl%nstep ! restart set this to nmax + tl%nstep
          tl%nstep=0
          tl%nstep0=2   
          ! copy prognostic variables:  tl%n0 into tl%nm1
          do ie=nets,nete
             elem(ie)%state%v(:,:,:,:,tl%nm1)=elem(ie)%state%v(:,:,:,:,tl%n0)
             elem(ie)%state%T(:,:,:,tl%nm1)=elem(ie)%state%T(:,:,:,tl%n0)
             elem(ie)%state%ps_v(:,:,tl%nm1)=elem(ie)%state%ps_v(:,:,tl%n0)
             elem(ie)%state%lnps(:,:,tl%nm1)=elem(ie)%state%lnps(:,:,tl%n0)
             elem(ie)%state%Q(:,:,:,:,tl%nm1)=elem(ie)%state%Q(:,:,:,:,tl%n0)
          enddo
       endif ! runtype==2

       if (hybrid%masterthread) then 
          write(iulog,*) "initial state from restart file:"
          write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
       end if
       
    else  ! initial run  RUNTYPE=0
       ! ===========================================================
       ! Initial Run  - compute initial condition
       ! ===========================================================
       if (hybrid%masterthread) then
          write(iulog,*) ' runtype: INITIAL primitive equations'
       endif
       ! ========================================================
       ! Initialize the test cases
       ! ========================================================
       if (test_case(1:10) == "baroclinic") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing Polvani-Scott-Thomas baroclinic instability test'
          end if
          
          call binst_init_state(elem, hybrid,nets,nete,hvcoord)
       else if (test_case(1:16) == "asp_gravity_wave") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing ASP gravity wave test'
          end if
          call asp_gravity_wave(elem, hybrid,hvcoord,nets,nete, sub_case)
       else if (test_case(1:12) == "asp_mountain") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing ASP mountain Rossby test'
          end if
          call asp_mountain(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:10) == "asp_rossby") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing ASP Rossby Haurwitz test'
          end if
          call asp_rossby(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:10) == "asp_tracer") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing pure tracer advection tests'
          end if
          call asp_tracer(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:14) == "asp_baroclinic") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing Jablonowski and Williamson ASP baroclinic instability test'
          end if
          call asp_baroclinic(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:13) == "jw_baroclinic") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing Jablonowski and Williamson baroclinic instability test V1'
          end if
          call jw_baroclinic(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:12) == "held_suarez0") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing Held-Suarez primitive equations test'
          end if
          call hs0_init_state(elem, hvcoord,nets,nete,300.0_real_kind)
       else if (test_case(1:10) == "aquaplanet") then
          if (hybrid%masterthread) then
             write(iulog,*)  'Initializing aqua planet with MJO-type forcing'
          end if
          if(moisture.eq."dry") then
             call binst_init_state(elem, hybrid,nets,nete,hvcoord)
          end if
          call aquaplanet_init_state(elem, hybrid,hvcoord,nets,nete,integration)
       else
          call abortmp('Error: unrecognized test case') 
       endif
       
       if (hybrid%masterthread) then
          write(iulog,*) '...done'
       end if
       
       ! scale PS to achieve prescribed dry mass
       call prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)  
       
       ! ========================================
       ! Print state and movie output
       ! ========================================
       
       if (hybrid%masterthread) then 
          write(iulog,*) "initial state:"
          write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
       end if
    end if  ! runtype

#endif

    ! For new runs, and branch runs, convert state variable to (Qdp)
    ! because initial conditon reads in Q, not Qdp
    ! restart runs will read dpQ from restart file
    ! need to check what CAM does on a branch run
    if (runtype==0 .or. runtype==2) then
       do ie=nets,nete
          elem(ie)%derived%omega_p(:,:,:) = 0D0
       end do

       if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
          do ie=nets,nete
             do t=1,3
             do q=1,qsize       
             do k=1,nlev
             do i=1,nv
             do j=1,nv          
                dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                     ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,t)
                elem(ie)%state%Qdp(i,j,k,q,t)=elem(ie)%state%Q(i,j,k,q,t)*dp  
             enddo
             enddo
             enddo
             enddo
             enddo
          enddo
       endif
    endif

#ifdef CAM
    if (hybrid%masterthread) then 
       ! CAM has set tstep based on dtime before calling prim_init2(), 
       ! so only now does HOMME learn the timstep.  print them out:
       if (phys_tscale/=0) then
          write(iulog,'(a,2f9.2)') "CAM physics timescale:       ",phys_tscale
       else
          write(iulog,'(a,2f9.2)') "CAM physics timescale: dtime"
       endif
       write(iulog,'(a,2f9.2)') "CAM dtime (dt_phys):         ",tstep*homme_nsplit*qsplit
       write(iulog,'(a,2f9.2)') "CAM dt_tracer, per RK stage: ",tstep*qsplit,(tstep*qsplit)/(rk_stage_user-1)
       write(iulog,'(a,2f9.2)') "CAM dt_dyn, per RK stage:    ",tstep*qsplit,tstep
 

       write(iulog,*) "initial state from CAM:"
       write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
    end if
#endif
    call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)

  end subroutine prim_init2




  subroutine leapfrog_bootstrap(elem,hybrid,nets,nete,tstep,tl,hvcoord)
  !
  ! leapfrog bootstrap code.  
  !
  ! take the equivilent of one timestep, but do it with a 
  ! dt/2 euler and a dt/2 leapfrog step  
  !
  use element_mod, only : element_t
  use hybvcoord_mod, only : hvcoord_t
  use time_mod, only : TimeLevel_t, time_at, timelevel_update
  use dimensions_mod, only :  nv, np, nlev, ne, qsize
  use control_mod, only : tstep_type

  type (element_t) , intent(inout)        :: elem(:)
  type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
  type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct
  integer, intent(in)                     :: nets  ! starting thread element number (private)
  integer, intent(in)                     :: nete  ! ending thread element number   (private)
  real(kind=real_kind), intent(in)        :: tstep          ! "timestep dependent" timestep
  type (TimeLevel_t), intent(inout)       :: tl


  ! local
  real(kind=real_kind) :: tstep_tmp,tstep_dyn
  integer :: i,ie


  tstep_dyn = tstep
  ! forward euler to get to tstep_dyn/2 (keep t=0 in nm1 timelevel) 
  ! (note: leapfrog tstep_dyn/4 with nm1=n0 is Euler with tstep_dyn/2 )
  tstep_tmp=tstep_dyn/4        

  call prim_run(elem, hybrid,nets,nete, tstep_tmp, tl, hvcoord, "forward")
  
  ! leapfrog with tstep_dyn/2 to get to tstep_dyn (keep t=0 in nm1 timelevel)
  tstep_tmp=tstep_dyn/2
  call prim_run(elem, hybrid,nets,nete, tstep_tmp, tl, hvcoord, "forward")  


  tl%nstep=tl%nstep-1        ! count all of that as 1 timestep

  end subroutine leapfrog_bootstrap





  subroutine prim_run(elem, hybrid,nets,nete, dt, tl, hvcoord, advance_name)
    use element_mod, only : element_t
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, time_at, timelevel_update, smooth
    use control_mod, only: statefreq, integration, tracer_advection_formulation,&
           TRACERADV_TOTAL_DIVERGENCE,TRACERADV_UGRADQ,ftype, energy_fixer, ftype, tstep_type
    use prim_advance_mod, only : prim_advance_exp, prim_advance_si, preq_robert3, &
         applycamforcing, applycamforcing_leapfrog
    use prim_advection_mod, only : Prim_Advec_Tracers_remap_rk2, prim_advec_tracers_lf
    use prim_state_mod, only : prim_printstate, prim_diag_scalars, prim_energy_halftimes
    use dimensions_mod, only : qsize,nlev,nv
    use parallel_mod, only : abortmp
#ifndef CAM
    use column_model_mod, only : ApplyColumnModel
#endif

    type (element_t) , intent(inout)        :: elem(:)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=real_kind), intent(in)        :: dt              ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    character(len=*), intent(in) :: advance_name
    real(kind=real_kind) :: st, st1, dp
    integer :: ie, t, q,k,i,j


    logical :: compute_diagnostics, compute_energy

    ! ===================================
    ! Main timestepping loop
    ! ===================================

    ! compute diagnostics and energy for STDOUT 
    ! compute energy if we are using an energy fixer

    compute_diagnostics=.false.
    compute_energy=energy_fixer>0
    if (MODULO(tl%nstep+1,statefreq)==0 .or. tl%nstep+1==tl%nstep0) then
       compute_diagnostics=.true.  
       compute_energy = .true.
    endif
    tot_iter=0.0       


    ! Forcing options for testing CAM-HOMME energy balance:
    if (ftype == -1) then
       ! disable all forcing, but allow moisture:
       do ie=nets,nete
          elem(ie)%derived%FQ = 0
          elem(ie)%derived%FM = 0
          elem(ie)%derived%FT = 0
       enddo
    endif
    if (ftype == -2) then
       ! disable moisture, but allow dynamics forcing
       do ie=nets,nete
          elem(ie)%state%Q = 0
          elem(ie)%state%Qdp = 0
          elem(ie)%derived%FQ = 0
       enddo
    endif
    if (ftype == -3) then
       ! disable forcing & moisture
       do ie=nets,nete
          elem(ie)%state%Q = 0
          elem(ie)%state%Qdp = 0
          elem(ie)%derived%FQ = 0
          elem(ie)%derived%FM = 0
          elem(ie)%derived%FT = 0
       enddo
    endif

    ! =================================
    ! energy, dissipation rate diagnostics.  Uses data at t-1,t 
    ! to compute diagnostics at t - 0.5.  
    ! small error in the t+.5 terms because at this
    ! point only state variables at t-1 has been Robert filtered.  
    ! =================================
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,1,nets,nete)
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,1,nets,nete)

    ! ===============
    ! initialize mean flux accumulation variables
    ! save U(t) for use by tracers 
    ! ===============
    do ie=nets,nete
       elem(ie)%derived%eta_dot_dpdn=0
       elem(ie)%derived%omega_p=0
       elem(ie)%derived%vn0=elem(ie)%state%v(:,:,:,:,tl%n0)
    enddo

    ! ===============
    ! Dynamical Step  uses Q at tl%n0
    ! ===============
!$OMP BARRIER
    if (integration == "explicit") then
       call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
            flt , hybrid,             &
            dt, tl, nets, nete, compute_diagnostics,tl%n0)
    else if (integration == "semi_imp") then
       call prim_advance_si(elem, nets, nete, cg(hybrid%ithr), blkjac, red, &
            refstate, hvcoord, deriv(hybrid%ithr), flt, hybrid, tl, dt)
       tot_iter=tot_iter+cg(hybrid%ithr)%iter
    end if


    ! ===============
    ! Tracer Advection  Needs U,V at timelevel n0 and eta_dot_dpdn at timellevel n0
    ! and maybe timelevel np1 which was computed in dynamics step above.  
    ! ===============

    if (qsize>0) then
       call Prim_Advec_Tracers_lf(elem, deriv(hybrid%ithr),hvcoord,flt_advection,hybrid,&
            dt,tl,nets,nete,compute_diagnostics)
    endif

    ! =================================
    ! energy, dissipation rate diagnostics.  Uses data at t and t+1
    ! to compute diagnostics at t + 0.5.
    ! =================================
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,2,nets,nete)
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,2,nets,nete)

    ! ===================================
    ! Compute Forcing Tendencies from nm1 data (for PROCESS SPLIT)
    ! or np1 data (for TIMESPLIT) and add tendencies into soluiton at timelevel np1
    ! ===================================       
#ifdef CAM
    ! ftype==1 means forcing is applied in dp_coupling.F90
    if (ftype<=0) call ApplyCAMForcing_leapfrog(elem, hvcoord,tl%n0,tl%np1,dt,nets,nete)
#else
    call ApplyColumnModel(elem, hybrid,cm(hybrid%ithr),dt)
#endif
    ! measure the effects of forcing
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,3,nets,nete)
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,3,nets,nete)


    ! =================================
    ! timestep is complete.  Now apply robert filter to all prognostic variables
    ! =================================
    if (smooth/=0) &
       call preq_robert3(tl%nm1,tl%n0,tl%np1,elem,hvcoord,nets,nete)
    ! measure the effects of Robert filter
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,4,nets,nete)
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,4,nets,nete)


    if (energy_fixer > 0) then
       if ( .not. compute_energy) call abortmp("ERROR: energy fixer needs compute_energy=.true")
       if ( ftype==0) call abortmp("ERROR: energy fixer cannot be used with ftype=0")
       call prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete)
       ! recompute fixed energy, if we are printing diagnostics:
       if (compute_diagnostics) call prim_energy_halftimes(elem,hvcoord,tl,4,nets,nete)
    endif

    ! =================================
    ! update dynamics time level pointers 
    ! =================================
    if(hybrid%ithr==0) call TimeLevel_update(tl,advance_name)

    ! ============================================================
    ! Print some diagnostic information 
    ! ============================================================
    !$OMP BARRIER
    if (compute_diagnostics) then
       if (hybrid%masterthread) then 
          write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
          if (integration == "semi_imp") write(iulog,*) "cg its=",cg(0)%iter
       end if
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    end if
  end subroutine prim_run




  subroutine prim_run_subcycle(elem, hybrid,nets,nete, dt, tl, hvcoord)
    use element_mod, only : element_t
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, time_at, timelevel_update, smooth, ptimelevels
    use control_mod, only: statefreq, integration, tracer_advection_formulation,&
           TRACERADV_TOTAL_DIVERGENCE,TRACERADV_UGRADQ,ftype, energy_fixer, ftype, qsplit, &
           compute_mean_flux
    use prim_advance_mod, only : prim_advance_exp, prim_advance_si, preq_robert3, applycamforcing, &
                                 applycamforcing_dynamics, prim_advance_exp
    use prim_advection_mod, only : prim_advec_tracers_remap_rk2
    use prim_state_mod, only : prim_printstate, prim_diag_scalars, prim_energy_halftimes
    use dimensions_mod, only : qsize,nlev,nv
    use parallel_mod, only : abortmp

    type (element_t) , intent(inout)        :: elem(:)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=real_kind), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    real(kind=real_kind) :: st, st1, dp, dt_q
    integer :: ie, t, q,k,i,j,n, n_Q


    logical :: compute_diagnostics, compute_energy


    ! ===================================
    ! Main timestepping loop
    ! ===================================
    dt_q = dt*qsplit


    ! compute diagnostics and energy for STDOUT 
    ! compute energy if we are using an energy fixer
    compute_diagnostics=.false.
    compute_energy=energy_fixer > 0
    if (MODULO(tl%nstep+qsplit,statefreq)==0 .or. tl%nstep+qsplit==tl%nstep0) then
       compute_diagnostics=.true.  
       compute_energy = .true.
    endif


    if (ftype < 0) print *,'ERROR: subcyling ftype<0 not yet coded'

#ifdef CAM
    ! ftype=2  Q was adjusted by physics, but apply u,T forcing here
    ! ftype=1  forcing was applied time-split in CAM coupling layer
    ! ftype=0 means forcing apply here
    if (ftype==0) call ApplyCAMForcing(elem, hvcoord,tl%n0,dt_q,nets,nete)
    if (ftype==2) call ApplyCAMForcing_dynamics(elem, hvcoord,tl%n0,dt_q,nets,nete)
#endif

    ! ===============
    ! initialize mean flux accumulation variables and store dp at n0 for use by advection
    ! ===============
    do ie=nets,nete
      elem(ie)%derived%eta_dot_dpdn=0    ! not used by subcycling code
      if (compute_mean_flux==1) then
         elem(ie)%derived%vn0=0.0D0
      else
         elem(ie)%derived%vn0=elem(ie)%state%v(:,:,:,:,tl%n0)
      endif
      elem(ie)%derived%omega_p=0

      do k=1,nlev
	  do i=1,nv
	    do j=1,nv      
		elem(ie)%derived%dp(i,j,k)=( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
		      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,tl%n0)
	    enddo
	  enddo
      enddo

    enddo


    n_Q = tl%n0   ! use Q at timelevel n0 for all dynamics steps:

    ! E(1) Energy at start of timestep, diagnostics at t-dt/2  (using t-dt, t and Q(t))
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,1,nets,nete,n_Q)
    ! qmass and variance, using Q(n0),Qdp(n0)
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,1,nets,nete)

    ! ===============
    ! Dynamical Step 
    ! ===============

    call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
         flt , hybrid, dt, tl, nets, nete, compute_diagnostics,n_Q)

    do n=2,qsplit

       if(hybrid%ithr==0) call TimeLevel_update(tl,"leapfrog")
!$OMP BARRIER
       call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
            flt , hybrid, dt, tl, nets, nete, .false.,n_Q)
       ! defer final robert, & timelevel update until after Q update.

    enddo

    ! ===============
    ! Tracer Advection.  needs U(t) (saved above) and U(t+1), dp(t+1) (now in elem%state%np1 )
    ! ===============
    ! Prim_Advec_Tracers will advance from n0 -> np1, so make sure Q(n_Q) data is in the n0 timelevel:
    if ( n_Q /= tl%n0 ) then
      do ie=nets,nete
        elem(ie)%state%Q(:,:,:,:,tl%n0)  = elem(ie)%state%Q(:,:,:,:,n_Q)
        elem(ie)%state%Qdp(:,:,:,:,tl%n0)  = elem(ie)%state%Qdp(:,:,:,:,n_Q)
      enddo
    endif
    if (qsize>0) call Prim_Advec_Tracers_remap_rk2(elem, deriv(hybrid%ithr),hvcoord,flt_advection,hybrid,&
         dt_q,tl,nets,nete,compute_diagnostics)

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_q - 2*dt 
    !   u(n0)    dynamics at  t+dt_q - dt       (Not yet Robert-filtered)
    !   u(np1)   dynamics at  t+dt_q 
    !
    !   Q(nm1)   undefined
    !   Q(n0)    Q at t
    !   Q(np1)   Q at t+dt_q
    if (compute_diagnostics) then
       call prim_diag_scalars(elem,hvcoord,tl,2,nets,nete)
       call prim_energy_halftimes(elem,hvcoord,tl,2,nets,nete,tl%np1)
    endif



    if (compute_diagnostics) then
       ! qmass and variance, using Q(np1),Qdp(np1).  n=3 and 4 values are all the same
       call prim_diag_scalars(elem,hvcoord,tl,3,nets,nete)
       call prim_diag_scalars(elem,hvcoord,tl,4,nets,nete)
       ! E(3) = Energy at end of timestep, before Robert:   t + dt_q  - dt/2
       ! (computed using n0,np1 dynamics and Q(np1) )
       ! E(3) only used for stdout diagnostics
       call prim_energy_halftimes(elem,hvcoord,tl,3,nets,nete,tl%np1)
    endif

    if (smooth/=0) &
         call preq_robert3(tl%nm1,tl%n0,tl%np1,elem,hvcoord,nets,nete)
    ! E(4) = Energy at end of timestep, after Robert:   t + dt_q  - dt/2
    ! (computed using n0,np1 dynamics and Q(np1) ).  Needed for fixer and some diagnostics
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,4,nets,nete,tl%np1)

    if (energy_fixer > 0) then
       if ( .not. compute_energy) call abortmp("ERROR: energy fixer needs compute_energy=.true")
       call prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete)
       ! recompute fixed energy, if we are printing diagnostics:
       if (compute_diagnostics) call prim_energy_halftimes(elem,hvcoord,tl,4,nets,nete,tl%np1)
    endif

    ! =================================
    ! update dynamics time level pointers 
    ! =================================
    if(hybrid%ithr==0) call TimeLevel_update(tl,"leapfrog")

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_q - dt       (Robert-filtered)
    !   u(n0)    dynamics at  t+dt_q 
    !   u(np1)   undefined
    !
    !   Q(nm1)   Q at t
    !   Q(n0)    Q at t+dt_q
    !   Q(np1)   undefined
 

    ! ============================================================
    ! Print some diagnostic information 
    ! ============================================================
    !$OMP BARRIER
    if (compute_diagnostics) then
       if (hybrid%masterthread) then 
          write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
       end if
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    end if

  end subroutine prim_run_subcycle




  subroutine prim_finalize(hybrid)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    ! ==========================
    ! end of the hybrid program
    ! ==========================
  end subroutine prim_finalize





  subroutine prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete)
! 
! non-subcycle code:
!  Solution is given at times u(t-1),u(t),u(t+1)
!  E(n=1) = energy at time t-.5
!  E(n=2) = energy at time t+.5
!  E(n=3) = energy at time t+.5  after Forcing applied, if ftype=0
!  E(n=4) = energy at time t+.5  after Robert filter
! Subcycle code
!  Solution is given at times u(t-1),u(t),u(t+1)
!  E(n=1) = energy at time t
!  E(n=2) = energy at time t+1
!  E(n=3) = energy at time t+1
!  E(n=4) = energy at time t+1
!
!  Tnew = T(t+1) + beta
!
! After Robert filter:
! E4 = .5[ cp_star(t1)*dpt1*T(t2)  + cp_star(t2)*dpt2*T(t1)  ] + KE + PE
! E5 = .5[ cp_star(t1)*dpt1*(T(t2)+beta)  + cp_star(t2)*dpt2*T(t1)  ] + KE + PE
!
! E5 = E4 +  .5< cp_star(t1)*dpt1 > beta
!  We want:
!        E5 = E1 
!         
!        E4 + .5<dp cp_star > beta = E1 
!       .5<dp cp_star> beta = (E1-E4) 
!   
    use parallel_mod, only: global_shared_buf, global_shared_sum
    use kinds, only : real_kind
    use dimensions_mod, only : nv, np, nlev, nelemd
    use hybvcoord_mod, only : hvcoord_t
    use element_mod, only : element_t
    use physical_constants, only : Cp, cpwater_vapor,g,dd_pi
    use physics_mod, only : Virtual_Specific_Heat
    use time_mod, only : timelevel_t
    use control_mod, only : moisture, tracer_advection_formulation,traceradv_total_divergence,energy_fixer, ftype
    use hybvcoord_mod, only : hvcoord_t
#ifdef CAM
    use repro_sum_mod, only: repro_sum
#endif
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
    integer :: t1,t2,t3,n,nets,nete,t_beta
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(inout)       :: tl

    integer :: ie,k,i,j,nm_f
    real (kind=real_kind), dimension(nv,nv,nlev)  :: dp   ! delta pressure
    real (kind=real_kind), dimension(nv,nv)  :: suml,suml2,v1,v2
    real (kind=real_kind) :: cp_star,qval

    real (kind=real_kind) :: psum(nets:nete,5),psum_g(5),beta,scale
    logical :: use_cp_star

! energy_fixer
!     1          cp_star(t1)*dp(t1)*T(t2)         (staggered)
!     2          cp*dp(t1)*T(t2)                  (staggered)
!     3          cp_star(t2)*dp(t2)*T(t2)
!     4          cp*dp(t2)*T(t2)
!
    t1=tl%n0     ! timelevel for cp_star dp 
    t2=tl%np1    ! timelevel for T

    use_cp_star = (moisture /= "dry")
    if (energy_fixer==2) use_cp_star = .false.
    if (energy_fixer==4) use_cp_star = .false.
    
    t_beta = t1
    scale = 2.0
    if (energy_fixer==3 .or. energy_fixer==4) then
       t_beta=t2
       scale = 1.0
    endif
    ! with staggered-in-time formula:
    !    compute <dp cp_star> at t1   
    !    E(Tnew) =(< dp(t1) cp_star(t1) (T(t2)+beta) > + < dp(t2) cp_star(t2) T(t1) >)/2 
    !    E(Tnew)-E(t2) = < dp(t1) cp_star(t1) beta >/2
    ! with everyting at t2: 
    !    compute <dp cp_star> at t2   
    !    E(Tnew)-E(t2) = < dp(t2) cp_star(t2) beta >
    !

    psum = 0
    do ie=nets,nete

       do k=1,nlev
          dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,t_beta)
       enddo
       suml=0
       do k=1,nlev
          do i=1,nv
          do j=1,nv
             if(use_cp_star)  then
                if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
                   qval = elem(ie)%state%Qdp(i,j,k,1,t_beta)/dp(i,j,k)
                else
                   qval = elem(ie)%state%Q(i,j,k,1,t_beta)
                endif
                cp_star= Virtual_Specific_Heat(qval)
                suml(i,j) = suml(i,j) + cp_star*dp(i,j,k) 
             else
                suml(i,j) = suml(i,j) + cp*dp(i,j,k) 
             endif

          enddo
          enddo
       enddo
       psum(ie,5) = psum(ie,5) + SUM(suml(:,:)*elem(ie)%spheremv(:,:))
       do n=1,4
          if (use_cp_star ) then
             psum(ie,n) = psum(ie,n) + SUM(  elem(ie)%spheremv(:,:)*&
                                    (elem(ie)%accum%PEner(:,:,n) + &
                                    elem(ie)%accum%IEner(:,:,n) + &
                                    elem(ie)%accum%KEner(:,:,n) ) )
          else
             psum(ie,n) = psum(ie,n) + SUM(  elem(ie)%spheremv(:,:)*&
                                    (elem(ie)%accum%PEner(:,:,n) + &
              (elem(ie)%accum%IEner(:,:,n)-elem(ie)%accum%IEner_wet(:,:,n)) + &
                                    elem(ie)%accum%KEner(:,:,n) ) )
          endif
       enddo
    enddo
#ifdef CAM
    do ie=nets,nete
       do n=1,5
          global_shared_buf(ie,n) = psum(ie,n)
       enddo
    enddo
!$OMP BARRIER
!$OMP MASTER
    call repro_sum(global_shared_buf, global_shared_sum, nelemd, nelemd, 5, commid=hybrid%par%comm )
!$OMP END MASTER
!$OMP BARRIER
    do n=1,5
       psum_g(n) = global_shared_sum(n)
    enddo
#else
    do ie=nets+1,nete
       psum(nets,:) = psum(nets,:) + psum(ie,:)
    end do
    psum_g = parallelsum(psum(nets,:),5,hybrid)
#endif
    beta = scale*( psum_g(1)-psum_g(4) )/psum_g(5)
    do ie=nets,nete
       elem(ie)%state%T(:,:,:,t2) =  elem(ie)%state%T(:,:,:,t2) + beta
    enddo


    end subroutine prim_energy_fixer


    subroutine smooth_topo_datasets(phis,sghdyn,sgh30dyn,elem,hybrid,nets,nete)
    use dimensions_mod, only : nv, np, nlev
    use control_mod, only : smooth_phis_numcycle,smooth_sgh_numcycle
    use hybrid_mod, only : hybrid_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
    use bndry_mod, only : bndry_exchangev
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t , laplace_sphere_wk
    use viscosity_mod, only : biharmonic_wk
    use time_mod, only : TimeLevel_t
    use prim_advance_mod, only : smooth_phis
    implicit none
    
    real (kind=real_kind), intent(inout)   :: phis(nv,nv,nets:nete)
    real (kind=real_kind), intent(inout)   :: sghdyn(nv,nv,nets:nete)
    real (kind=real_kind), intent(inout)   :: sgh30dyn(nv,nv,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    integer , intent(in) :: nets,nete
    ! local
    integer :: ie
    real (kind=real_kind) :: minf 

    minf=-9e9
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to PHIS"
    call smooth_phis(phis,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_phis_numcycle)

    minf=0
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to SGH"
    call smooth_phis(sghdyn,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_sgh_numcycle)
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to SGH30"
    call smooth_phis(sgh30dyn,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_sgh_numcycle)

    end subroutine smooth_topo_datasets
    


end module prim_driver_mod



