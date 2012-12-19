#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module element_mod
  ! ------------------------------
  use kinds, only : real_kind, long_kind, int_kind
  ! ------------------------------
  use coordinate_systems_mod, only : spherical_polar_t, cartesian2D_t
  !--------------------------
  use dimensions_mod, only : np, nv, nvsq, nlev, nlevp, qsize_d
  !--------------------------
  use edge_mod, only : edgedescriptor_t, rotation_t
  !--------------------------
  use gridgraph_mod, only : gridvertex_t
  !--------------------------
  implicit none
  private
  integer, parameter :: timelevels=3

#ifdef _PRIM
  type, public :: elem_state_t
     sequence
!
!    note: variables (and sequence) must match that in prim_restart_mod.F90
!
! prognostic variables
!
     real (kind=real_kind) :: v(nv,nv,2,nlev,timelevels)   ! velocity                                            1
     real (kind=real_kind) :: T(np,np,nlev,timelevels)     ! temperature                                         2
     real (kind=real_kind) :: lnps(np,np,timelevels)       ! log surface pressure                                3
     real (kind=real_kind) :: ps_v(nv,nv,timelevels)       ! surface pressure on v grid                          4
     real (kind=real_kind) :: phis(np,np)         ! surface geopotential (prescribed)                            5
     ! qsize = 1 is related to the mixing ratio
     ! everything else are passive tracers that can eventually
     ! be forced by the column model.
     real (kind=real_kind) :: Q(nv,nv,nlev,qsize_d,timelevels)  ! Tracers             
     real (kind=real_kind) :: Qdp(nv,nv,nlev,qsize_d,timelevels)  ! Tracers           

  end type elem_state_t

  ! JPE: This parameter must match the number of variables in the state
  ! structure all of which are assumed to be of kind=real_kind.  This is a
  ! requirement for restart I/O.

  integer(kind=int_kind),public,parameter :: StateComponents=7


  type, public :: derived_state_t
     sequence


!    storage needed when subcycling tracers/dynamics 
     real (kind=real_kind) :: vn0(nv,nv,2,nlev)   ! velocity at n0 saved for use by tracers when sub-cycling
     real (kind=real_kind) :: vstar(nv,nv,2,nlev) ! velocity on Lagrangian surfaces

! diagnostic variables computed during explicit timestep
     real (kind=real_kind) :: phi(np,np,nlev)     ! geopotential                      
     real (kind=real_kind) :: omega_p(nv,nv,nlev) ! vertical tendency (derived)       
     real (kind=real_kind) :: eta_dot_dpdn(nv,nv,nlevp)  ! vertical velocity to be used for tracers
! Semi-Implicit diagnostic variables, computed during
! explict component, and then reused in Helmholtz component.
     real (kind=real_kind) :: grad_lnps(nv,nv,2)  ! gradient of log surface pressure               
     real (kind=real_kind) :: zeta(nv,nv,nlev)    ! relative vorticity                             
     real (kind=real_kind) :: div(nv,nv,nlev,timelevels)     ! divergence                          

! ---------------------- OG added
!these fields are for consistency and limiters
!field dp (same for tracers and dynamics at the beginning of each time step)
!is to construct dp_tracers at physics timestep
     real (kind=real_kind) :: dp(nv,nv,nlev)       
!divdp_proj is DSSed divdp
     real (kind=real_kind) :: divdp(nv,nv,nlev) 
     real (kind=real_kind) :: divdp_proj(nv,nv,nlev) 
!-------------------------------------------------------- 

     ! Primitive equations forcings
#ifdef CAM
     real (kind=real_kind) :: FQ(nv,nv,nlev,qsize_d, 1) ! F-Tracers   
     real (kind=real_kind) :: FM(nv,nv,2,nlev, 1)       ! F-momentum  
     real (kind=real_kind) :: FT(nv,nv,nlev, 1)         ! F-Temperature
#else
     ! when does forcing ever need more than 1 timelevel?  
     real (kind=real_kind) :: FQ(nv,nv,nlev,qsize_d, timelevels) ! F-Tracers  
     real (kind=real_kind) :: FM(nv,nv,2,nlev, timelevels) ! F-momentum       
     real (kind=real_kind) :: FT(nv,nv,nlev, timelevels)    ! F-Temperature   
#endif
     real (kind=real_kind) :: FQps(nv,nv,timelevels)         ! implied forcing of FQ on ps_v
                                                            ! used when conserving dry mass in the presence of precipitation 


  end type derived_state_t

  type, public :: elem_accum_t
     sequence
!     real (kind=real_kind) :: u(nv,nv,nlev)       ! zonal velocity on sphere
!     real (kind=real_kind) :: T(np,np,nlev)       ! temperature
!     real (kind=real_kind) :: Q(nv,nv,nlev) ! tracers
!     real (kind=real_kind) :: ke(nv,nv,nlev)      ! kinetic energy

!
! **** ENERGY DIAGNOSTICS ****    
!
! Energy equation:   
! KE_t  = T1 + T2  + D1   + Err   +  vertical & horizontal advection terms
! IE_t  = S1 + D2                 +  vertical & horizontal advection terms
! PE_t  = S2        
!
!  KEvert*  =  KE net vertical advection    (should be zero) 
!  KEhoriz* =  KE net horizonatl advection  (should be zero)
!  IEvert*  =  IE net vertical advection    (should be zero) 
!  IEhoriz* =  IE net horizonatl advection  (should be zero)
!
! With leapfrog, energy equations are all exact except KE 
! (has an Err term that goes to zero as dt**2)
!
! Transfer terms:
! T1   = -< dp/dn u, RT_v/p grad_p >     KE<->IE:   T1 + T2-T2_s = S1
! T2   = -< dp/dn u, grad_phi >          KE<->PE:   T2_s         = S2
! T2_s = -< dp/dn u, grad_phis >
! S1   = < Cp_star dp/dn , RT omega_p/Cp_star >  
! S2   = -< div (u dp/dn), phis >                
!
#ifdef ENERGY_DIAGNOSTICS
     real (kind=real_kind) :: KEvert1(nv,nv)      ! term from continuity equ
     real (kind=real_kind) :: KEvert2(nv,nv)      ! term from momentum equ
     real (kind=real_kind) :: IEvert1(nv,nv)      ! term from continuity equ
     real (kind=real_kind) :: IEvert2(nv,nv)      ! term from T equ
     real (kind=real_kind) :: IEvert1_wet(nv,nv)  ! wet term from continuity equ
     real (kind=real_kind) :: IEvert2_wet(nv,nv)  ! wet term from T equ

     real (kind=real_kind) :: KEhorz1(nv,nv)      ! at time t
     real (kind=real_kind) :: KEhorz2(nv,nv)      ! after calling time_advance, these will be at time t-1
     real (kind=real_kind) :: IEhorz1(nv,nv)    
     real (kind=real_kind) :: IEhorz2(nv,nv)    
     real (kind=real_kind) :: IEhorz1_wet(nv,nv)    
     real (kind=real_kind) :: IEhorz2_wet(nv,nv)    

     real (kind=real_kind) :: T1(nv,nv)    
     real (kind=real_kind) :: T2(nv,nv) 
     real (kind=real_kind) :: T2_s(nv,nv)    
     real (kind=real_kind) :: S1(nv,nv)    
     real (kind=real_kind) :: S1_wet(nv,nv)    
     real (kind=real_kind) :: S2(nv,nv)    

     ! the KE conversion term and diffusion term
     real (kind=real_kind) :: DIFF(nv,nv,2,nlev) ! net hypervis term
     real (kind=real_kind) :: DIFFT(nv,nv,nlev) ! net hypervis term
     real (kind=real_kind) :: CONV(nv,nv,2,nlev) ! dpdn u dot CONV = T1 + T2
#endif
     ! the "4" timelevels represents data computed at:
     !  1  t-.5   
     !  2  t+.5   after dynamics
     !  3  t+.5   after forcing
     !  4  t+.5   after Robert
     ! after calling TimeLevelUpdate, all time above decrease by 1.0
     real (kind=real_kind) :: KEner(nv,nv,4)       
     real (kind=real_kind) :: PEner(nv,nv,4)       
     real (kind=real_kind) :: IEner(nv,nv,4)    
     real (kind=real_kind) :: IEner_wet(nv,nv,4)    
     real (kind=real_kind) :: Qvar(nv,nv,qsize_d,4)  ! Q variance at half time levels   
     real (kind=real_kind) :: Qmass(nv,nv,qsize_d,4) ! Q mass at half time levels
     real (kind=real_kind) :: Q1mass(nv,nv,qsize_d)  ! Q mass at full time levels
     real (kind=real_kind) :: mass_added(qsize_d)    ! mass added by qneg fixer

  end type elem_accum_t

#elif defined _PRIMDG
  ! definition of elem_state_t for DG Primitive Equations version of model
  type, public :: elem_state_t
     sequence
     real (kind=real_kind) :: p(np,np,nlev,timelevels)
     real (kind=real_kind) :: phis(np,np)                  ! surface geopotential
     real (kind=real_kind) :: gradps(nv,nv,2)              ! gradient of surface geopotential
     real (kind=real_kind) :: v(nv,nv,2,nlev,timelevels)   ! contravarient comp
     !=======================================================================================================!
     real (kind=real_kind) :: couv(nv,nv,2,nlev) !covariant velocities
     real (kind=real_kind) :: uv(nv,nv,2,nlev)   !redundant copy of v (eliminate?)
     real (kind=real_kind) :: uv0(nv,nv,2,nlev)  !temp variable for velocities (eliminate?)
     real (kind=real_kind) :: pgrads(nv,nv,2,nlev)
     !=======================================================================================================!
     real (kind=real_kind) :: psi(nv,nv,nlev)
     real (kind=real_kind) :: phi(nv,nv,nlev)
     real (kind=real_kind) :: ht(nv,nv,nlev)
     real (kind=real_kind) :: T(nv,nv,nlev,timelevels)
     real (kind=real_kind) :: pt3d(nv,nv,nlev)
     real (kind=real_kind) :: peta(nv,nv,nlev)
     real (kind=real_kind) :: dp3d(nv,nv,nlev)
     !=======================================================================================================!
     real (kind=real_kind) :: zeta(nv,nv,nlev)
     real (kind=real_kind) :: pr3d(nv,nv,nlev+1)
     real (kind=real_kind) :: pr3d_ref(nv,nv,nlev+1)
     real (kind=real_kind) :: gp3d(nv,nv,nlev+1)
     !=======================================================================================================!
     real (kind=real_kind) :: ptop(nv,nv)
     real (kind=real_kind) :: sgp(nv,nv)
     real (kind=real_kind) :: tbar(nlev)
     !=======================================================================================================!
  end type elem_state_t

  ! this type is just a place holder at the moment.  In the future some of the fields in elem_state_t
  ! should be moved into derived_state_t to mirror the _PRIM verion
  type, public :: derived_state_t
     sequence
     real (kind=real_kind) :: dummmy
  end type derived_state_t


  type, public :: elem_accum_t
    sequence
    real (kind=real_kind) :: u(nv,nv,nlev)       ! zonal velocity on sphere
    real (kind=real_kind) :: T(nv,nv,nlev)       ! temperature
    real (kind=real_kind) :: ke(nv,nv,nlev)      ! kinetic energy
  end type

#else
  !
  ! SHALLOW WATER STRUCT
  !
  type, public :: elem_state_t
     sequence
     real (kind=real_kind) :: p(np,np,nlev,timelevels)
     real (kind=real_kind) :: ps(np,np)           ! surface geopotential
     real (kind=real_kind) :: gradps(nv,nv,2)     ! gradient of surface geopotential     
     real (kind=real_kind) :: v(nv,nv,2,nlev,timelevels)   ! contravarient comp
     !=======================================================================================================!
     ! This are all fields added for DG (Shallow-Water)							!
     !=======================================================================================================!
#ifdef _SWDG     
     real (kind=real_kind) :: couv(nv,nv,2,nlev) 	! couvarient velocity
     real (kind=real_kind) :: psi(nv,nv,nlev)     	! (ht+hs)*metdet(G)
     real (kind=real_kind) :: ht(nv,nv,nlev)      	! height variable     
     real (kind=real_kind) :: hs(nv,nv,nlev) 	 	! mountain height
#endif
     !=======================================================================================================!
  end type elem_state_t

  type, public :: derived_state_t
     sequence
     !=======================================================================================================!
     ! all shallow water diagnostic variables should be moved here.
     ! %state is for prognostic variables
     !=======================================================================================================!
     real (kind=real_kind) :: dummy
  end type derived_state_t

!=======================================================================================================! 
#endif

  type, public :: index_t
     sequence
     integer(kind=int_kind) :: ia(nvsq),ja(nvsq)
     integer(kind=int_kind) :: is,ie
     integer(kind=int_kind) :: NumUniquePts
     integer(kind=int_kind) :: UniquePtOffset
  end type index_t

  type, public :: element_t
     sequence

     integer(kind=int_kind) :: LocalId
     integer(kind=int_kind) :: GlobalId
     !=====================================
     !Add the link list hooks
     !=====================================
     type (element_t), pointer :: prev
     type (element_t), pointer :: next

     ! Coordinate values of element points

     type (spherical_polar_t) :: spherev(nv,nv)     ! Spherical coordinates of
     ! element points on velocity grid
     type (spherical_polar_t) :: spherep(np,np)     ! Spherical coordinates of
     ! element points on pressure grid
     type (cartesian2D_t)     :: cartv(nv,nv)       ! Cartesian coordinates of 
     ! element points on velocity points
     type (cartesian2D_t)     :: cartp(np,np)       ! Cartesian coordinates of
     ! element points on pressure points

     real (kind=real_kind)    :: dx                 ! Size of element x dimension   
     ! (assumes rectilinear)
     real (kind=real_kind)    :: dy                 ! Size of element in y direction 
     ! (assumes rectilinear)

     ! Edge connectivity information

     type (GridVertex_t)       :: vertex                 ! Element grid vertex information
     type (EdgeDescriptor_t) :: desc

     type (elem_state_t)       :: state

     type (derived_state_t)    :: derived
#if defined _PRIM || defined _PRIMDG
     type (elem_accum_t)       :: accum
#endif
     ! Metric terms 

     real (kind=real_kind)    :: met(2,2,nv,nv)      ! metric tensor on velocity grid
     real (kind=real_kind)    :: metinv(2,2,nv,nv)   ! metric tensor on velocity grid
     real (kind=real_kind)    :: metdet(nv,nv)       ! g = SQRT(det(g_ij)) on velocity
     real (kind=real_kind)    :: metdetp(np,np)      ! g = SQRT(det(g_ij)) on pressure grid
     real (kind=real_kind)    :: rmetdetp(np,np)     ! 1/metdet on pressure grid
     real (kind=real_kind)    :: D(2,2,nv,nv)        ! Map covariant field on cube to vector field on the sphere
     real (kind=real_kind)    :: Dinv(2,2,nv,nv)     ! Map vector field on the sphere to covariant v on cube

     ! Mass matrix terms for an element on a cube face

     real (kind=real_kind)    :: mp(np,np)          ! mass matrix on pressure grid
     real (kind=real_kind)    :: mv(nv,nv)          ! mass matrix on velocity grid
     real (kind=real_kind)    :: rmv(nv,nv)         ! inverse mass matrix on velocity grid

     ! Mass matrix terms for an element on the sphere
     ! This mass matrix is used when solving the equations in weak form
     ! with the natural (surface area of the sphere) inner product
     real (kind=real_kind)    :: spheremv(nv,nv)          ! mass matrix on velocity grid
     real (kind=real_kind)    :: rspheremv(nv,nv)         ! inverse mass matrix on velocity grid

     integer(kind=long_kind)  :: gdofV(nv,nv)        ! Global degree of freedom (V-grid) 
     integer(kind=long_kind)  :: gdofP(np,np)        ! Global degree of freedom (P-grid)

     ! Coreolis term

     real (kind=real_kind)    :: fcor(nv,nv)        ! coreolis term

     ! Solver weights (used only for non-staggered grid

     real (kind=real_kind)    :: solver_wts(nv,nv)

     type (index_t) :: idxP
     type (index_t),pointer :: idxV
     integer :: FaceNum

     ! force element_t to be a multiple of 8 bytes.  
     ! on BGP, code will crash (signal 7, or signal 15) if 8 byte alignment is off
     ! check core file for:
     ! core.63:Generated by interrupt..(Alignment Exception DEAR=0xa1ef671c ESR=0x01800000 CCR0=0x4800a002)
     integer :: dummy
  end type element_t


  type, public :: eroot_t
     sequence
     type(element_t), pointer :: first
  end type eroot_t

  integer, public :: NumEdges

  type (eroot_t), public :: eroot

  public :: element_coordinates

  public :: LLAddEdge,LLFindEdge, LLInsertEdge
  public :: LLSetEdgeCount,LLGetEdgeCount
  public :: LLFree


  public :: GetColumnIdP,GetColumnIdV
#if 0
  interface assignment  ( = ) 
     module procedure copy_node
  end interface
#endif

contains

  ! =======================================
  !  GetColumnIdP:
  !  
  !  Gives a unique identifier for a Physics 
  !  column on the P-grid
  ! =======================================
  function GetColumnIdP(elem,i,j) result(col_id)
    type(element_t), intent(in) :: elem
    integer, intent(in) :: i,j
    integer :: col_id

    col_id = elem%gdofP(i,j)

  end function GetColumnIdP
  ! =======================================
  !  GetColumnIdV:
  !  
  !  Gives a unique identifier for a Physics 
  !  column on the V-grid
  ! =======================================
  function GetColumnIdV(elem,i,j) result(col_id)
    type(element_t), intent(in) :: elem
    integer, intent(in) :: i,j
    integer :: col_id

    col_id = elem%gdofV(i,j)

  end function GetColumnIdV

  ! =======================================
  ! element_coordinates:
  !
  ! Initialize 2D rectilinear element 
  ! colocation points
  !
  ! =======================================

  function element_coordinates(start,end,points) result(cart)
    use kinds, only : longdouble_kind
    type (cartesian2D_t), intent(in) :: start
    type (cartesian2D_t), intent(in) :: end
    real (kind=longdouble_kind), intent(in) :: points(:)

    type (cartesian2D_t) :: cart(SIZE(points),SIZE(points))

    ! Local variables

    integer i,j

    do j=1,SIZE(points)
       do i=1,SIZE(points)
          cart(i,j)%x=0.50D0*((start%x+end%x)+(end%x-start%x)*points(i))
          cart(i,j)%y=0.50D0*((start%y+end%y)+(end%y-start%y)*points(j))
       end do
    end do

  end function element_coordinates

  subroutine LLSetEdgeCount(value)
    implicit none
    integer,intent(in)   :: value
    NumEdges=value
  end subroutine LLSetEdgeCount

  subroutine LLGetEdgeCount(value)
    implicit none
    integer,intent(out)  :: value
    value=NumEdges
  end subroutine LLGetEdgeCount

  recursive subroutine copy_node(node2,node1)

    type (element_t), intent(out) :: node2
    type (element_t), intent(in)  :: node1

    node2%LocalId = node1%LocalId
    node2%GlobalId = node1%GlobalId
    node2%prev       = node1%prev
    node2%next       = node1%next

  end subroutine copy_node

  subroutine LLFree(List)

    implicit none
    type(eroot_t) :: List
    type(element_t), pointer :: temp_node
    integer :: nlist,i


    temp_node => List%first
    ! Find the end of the list
    do while(associated(temp_node%next))
       temp_node => temp_node%next
    enddo

    temp_node => temp_node%prev
    !Now step back and deallocate all entries  
    do while(associated(temp_node))
       deallocate(temp_node%next)
       temp_node => temp_node%prev
    enddo

  end subroutine LLFree

  subroutine LLInsertEdge(EdgeList,Gid,Lid)
    type (eroot_t), intent(inout) :: EdgeList
    integer, intent(in)  :: Gid
    integer, intent(inout) :: Lid
    logical :: found

    call LLFindEdge(EdgeList,Gid,found) 
    if(.not. found) then 
       call LLAddEdge(EdgeList,Gid,Lid) 
    endif

  end subroutine LLInsertEdge

  subroutine LLFindEdge(Edge,Gid,found)

    type (eroot_t), intent(in) :: Edge
    integer, intent(in)  :: Gid
    logical, intent(out) :: found

    type (element_t), pointer :: temp_node

    found =.FALSE.

    temp_node => Edge%first
    do while(associated(temp_node) .and. (.not. found))
       if(Gid .eq. temp_node%GlobalId) then 
          found = .TRUE. 
       else
          temp_node => temp_node%next
       endif
    enddo
  end subroutine LLFindEdge

  subroutine LLAddEdge(EdgeList,Gid,Lid)
    type (eroot_t), intent(inout) :: EdgeList
    integer, intent(in)  :: Gid
    integer, intent(out)  :: Lid

    type(element_t), pointer :: temp_node
    type(element_t), pointer  :: new_node
    type(element_t), pointer :: parent

    temp_node => EdgeList%first
    parent    => EdgeList%first

    do while(associated(temp_node))
       parent => temp_node
       temp_node => parent%next
    enddo
    allocate(new_node)
    NumEdges = NumEdges + 1

    new_node%GlobalId=Gid
    new_node%LocalId=NumEdges
    NULLIFY(new_node%next)
    new_node%prev => parent

    if(associated(EdgeList%first)) then
       parent%next => new_node 
    else
       EdgeList%first => new_node 
    endif
    Lid = NumEdges

  end subroutine LLAddEdge

end module element_mod
