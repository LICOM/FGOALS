#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module GridGraph_mod
  !-------------------------
  use kinds, only : real_kind, iulog
  !-------------------------
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  !-----
  implicit none


  private


  type, public :: GridVertex_t
      sequence
      integer                   :: nbrs(8)   
      integer                   :: wgtV(8)   ! The weights for edges defined by neighbors array
      integer                   :: wgtP(8)   ! The weights for edges defined by neighbors array
      integer                   :: degree   ! degree of grid element (number of neighbors)
      integer                   :: number   ! element number
      integer                   :: partition 
      integer                   :: SpaceCurve  ! index in Space-Filling curve
  end type GridVertex_t

  type, public :: EdgeIndex_t
      sequence
      integer, pointer            :: ixV(:)
      integer, pointer            :: ixP(:)
      integer, pointer            :: iyV(:)
      integer, pointer            :: iyP(:)
  end type EdgeIndex_t

  type, public :: GridEdge_t
      sequence
      logical                      :: reverse
      integer                      :: head_face  ! needed if head vertex has shape (i.e. square)
      integer                      :: tail_face  ! needed if tail vertex has shape (i.e. square)
      type (GridVertex_t),pointer  :: head       ! edge head vertex
      type (GridVertex_t),pointer  :: tail       ! edge tail vertex

#ifdef TESTGRID
      integer                      :: wgtV        ! amount of information which must be transfered
      integer                      :: wgtP        ! amount of information which must be transfered
      type (EdgeIndex_t)           :: HeadIndex
      type (EdgeIndex_t)           :: TailIndex
#endif
  end type GridEdge_t
  
! ==========================================
! Public Interfaces
! ==========================================

  public :: set_GridVertex_number
  public :: PrintGridVertex
 
  public :: initgridedge
  public :: gridedge_search
  public :: gridedge_type
  public :: grid_edge_uses_vertex
  public :: PrintGridEdge
  public :: CheckGridNeighbors
  public :: PrintChecksum

  public :: CreateSubGridGraph
  public :: FreeGraph

  interface assignment ( = )
      module procedure copy_gridedge
      module procedure copy_edgeindex
      module procedure copy_gridvertex
  end interface

contains

! =====================================
! copy edge:
! copy device for overloading = sign.
! =====================================

  recursive subroutine copy_gridedge(edge2,edge1)

    type (GridEdge_t), intent(out) :: edge2
    type (GridEdge_t), intent(in)  :: edge1


    edge2%tail_face = edge1%tail_face
    edge2%head_face = edge1%head_face
    edge2%reverse   = edge1%reverse

    if (associated(edge1%tail)) then
       edge2%tail=>edge1%tail
    end if
    if (associated(edge1%head)) then
       edge2%head=>edge1%head
    end if

#ifdef TESTGRID
    edge2%wgtV       = edge1%wgtV
    edge2%wgtP       = edge1%wgtP

    edge2%TailIndex = edge1%TailIndex
    edge2%HeadIndex = edge1%HeadIndex
#endif

  end subroutine copy_gridedge

  recursive subroutine copy_gridvertex(vertex2,vertex1)
	
    implicit none 

    type (GridVertex_t), intent(out)   :: vertex2
    type (GridVertex_t), intent(in)    :: vertex1

    integer                            :: i,n
   
    
!JMD     vertex2%size      = vertex1%size
 
     n = SIZE(vertex1%nbrs)
     do i=1,n
	vertex2%nbrs(i) = vertex1%nbrs(i)
        vertex2%wgtV(i)  = vertex1%wgtV(i)
        vertex2%wgtP(i)  = vertex1%wgtP(i)
     enddo

     vertex2%degree     = vertex1%degree
     vertex2%number     = vertex1%number
     vertex2%partition  = vertex1%partition 
     vertex2%SpaceCurve = vertex1%SpaceCurve 

  end subroutine copy_gridvertex

  recursive subroutine copy_edgeindex(index2,index1)
  
  type (EdgeIndex_t), intent(out) :: index2
  type (EdgeIndex_t), intent(in)  :: index1

  if(associated(index1%iyV)) then 
    index2%iyV => index1%iyV
  endif

  if(associated(index1%iyP)) then 
    index2%iyP => index1%iyP
  endif

  if(associated(index1%ixV)) then 
     index2%ixV => index1%ixV
  endif

  if(associated(index1%ixP)) then 
     index2%ixP => index1%ixP
  endif
 
  end subroutine copy_edgeindex

  subroutine FreeGraph(Vertex)

     implicit none
     type (GridVertex_t)           :: Vertex(:)
     integer                       :: i,nelem

     nelem = SIZE(Vertex)

!JMD     do i=1,nelem
!JMD        deallocate(Vertex(i)%wgtV)
!JMD        deallocate(Vertex(i)%wgtG)
!JMD        deallocate(Vertex(i)%nbrs)
!JMD     enddo

  end subroutine FreeGraph

!===========================
! search edge list for match
!===========================

  function gridedge_search(nvert1,nvert2,edge) result(number)

    integer, intent(in) :: nvert1
    integer, intent(in) :: nvert2
    type(GridEdge_t), intent(in) :: edge(:)
    integer :: number

    integer :: tmp
    integer :: head
    integer :: tail

    integer :: nedge
    integer :: i

    nedge=SIZE(edge)

    tail=nvert1
    head=nvert2

    if (tail > head) then
       tmp  = tail
       tail = head
       head = tmp
    end if

    do i=1,nedge
       if (edge(i)%tail%number==tail .and. edge(i)%head%number==head)then
          number=i
       end if
    end do

  end function gridedge_search


  function gridedge_type(edge) result(type)

    use params_mod, only : INTERNAL_EDGE, EXTERNAL_EDGE
    type (GridEdge_t), intent(in)  :: edge
    integer                        :: type

    if (edge%head%partition==edge%tail%partition) then
        type=INTERNAL_EDGE
    else
        type=EXTERNAL_EDGE
    endif

  end function gridedge_type







  function grid_edge_uses_vertex(Vertex,Edge) result(log)

    type(GridVertex_t), intent(in) :: Vertex
    type(GridEdge_t),   intent(in) :: Edge
    logical :: log
    integer  :: number

    number = Vertex%number
    if(number == Edge%head%number .or. number == Edge%tail%number) then
        log = .TRUE.
    else
        log = .FALSE.
    endif

  end function grid_edge_uses_vertex

  subroutine PrintChecksum(TestPattern,Checksum)

   use dimensions_mod, only : nlev, nelemd, nv

   implicit none

   real(kind=real_kind), target,intent(in)   :: TestPattern(:,:,:,:)
   real(kind=real_kind), target,intent(in)   :: Checksum(:,:,:,:)

   integer                                  :: i,k,ix,iy

   print *
   write (iulog,*) 'checksums:'
   do i=1,nelemd
     !  Lets start out only looking at the first element
        write(iulog,*)
        do k=1,nlev
        do iy=1,nv
        do ix=1,nv
           write(iulog,*)INT(TestPattern(ix,iy,k,i))," checksum = ",INT(Checksum(ix,iy,k,i))
        enddo
        enddo
        enddo
   enddo


  end subroutine PrintChecksum

  subroutine CreateSubGridGraph(Vertex,SVertex,local2global)

    implicit none
	
    type (GridVertex_t),intent(in)         :: Vertex(:)
    type (GridVertex_t),intent(inout)      :: SVertex(:)
    integer,intent(in)                     :: local2global(:)

    integer                                :: nelem,nelem_s,n,ncount
    integer                                :: inbr,i,ig,j
    
    integer,allocatable                    :: global2local(:)
    logical, parameter    :: Debug = .FALSE.


    nelem   = SIZE(Vertex)
    nelem_s = SiZE(SVertex) 

    if(Debug) write(iulog,*)'CreateSubGridGraph: point #1'
    allocate(global2local(nelem))
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #2'

    global2local(:) = 0
    do i=1,nelem_s
        ig = local2global(i)
	global2local(ig) = i
    enddo
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #3'
    do i=1,nelem_s
	ig = local2global(i)

        if(Debug) write(iulog,*)'CreateSubGridGraph: point #4'
	call copy_gridvertex(SVertex(i),Vertex(ig))
	n = SIZE(SVertex(i)%nbrs(:))
	! ==============================================
        ! Apply the correction to the neighbors list to 
        ! reflect new subgraph numbers 
	! ==============================================
        ncount=0
        if(Debug) write(iulog,*)'CreateSubGridGraph: point #5'
        do j=1,n
           if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.1 size(global2local) Svertex(i)%nbrs(j) ', &
			size(global2local), Svertex(i)%nbrs(j)
           if(Svertex(i)%nbrs(j) > 0) then 
	      inbr = global2local(Svertex(i)%nbrs(j))
              if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.2'
	      if(inbr .gt. 0) then 
	         Svertex(i)%nbrs(j) = inbr
	         ncount = ncount+1
	      else 
	         Svertex(i)%wgtV(j)  = 0	
	         Svertex(i)%wgtP(j)  = 0	
	         Svertex(i)%nbrs(j) = -1
	      endif
           endif
           if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.3'
        enddo
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #6'
        Svertex(i)%number = i
     enddo
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #7'
     deallocate(global2local)
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #8'

  end subroutine CreateSubGridGraph

  subroutine PrintGridEdge(Edge)

    implicit none
    type (GridEdge_t), intent(in) :: Edge(:)

    integer           :: i,nedge,ii,wgtV

    nedge = SIZE(Edge)

    write(iulog,95)
    do i=1,nedge
          ii=Edge(i)%tail_face
          wgtV=Edge(i)%tail%wgtV(ii)       
          write(iulog,100) i, &
               Edge(i)%tail%number,Edge(i)%tail_face, wgtV, &
               Edge(i)%head%number,Edge(i)%head_face, gridedge_type(Edge(i))
    enddo
  95 format(5x,'GRIDEDGE #',3x,'Tail (face)',5x,'Head (face)',3x,'Type')
 100 format(10x,I6,8x,I4,1x,'(',I1,')  --',I2,'--> ',I6,1x,'(',I1,')',5x,'[',I1,']')

  end subroutine PrintGridEdge

! ==========================================
! set_GridVertex_neighbors:
!
! Set global element number for element elem
! ==========================================

  subroutine set_GridVertex_number(elem,number)

    type(GridVertex_t)         :: elem
    integer                 :: number

    elem%number=number

  end subroutine set_GridVertex_number

  subroutine PrintGridVertex(Vertex)

    implicit none 
    type (GridVertex_t), intent(in),target :: Vertex(:)
  
    integer        :: i,nvert

    nvert = SIZE(Vertex)
	
    write(iulog,98)
    do i=1,nvert
#if 0
    write(iulog,*) Vertex(i)%number, Vertex(i)%partition
#else
       write(iulog,99) Vertex(i)%number, Vertex(i)%partition, Vertex(i)%degree, &
	   Vertex(i)%nbrs(west),Vertex(i)%wgtV(west), &
           Vertex(i)%nbrs(east),Vertex(i)%wgtV(east), &
           Vertex(i)%nbrs(south),Vertex(i)%wgtV(south), &
           Vertex(i)%nbrs(north),Vertex(i)%wgtV(north), &
           Vertex(i)%nbrs(swest),Vertex(i)%wgtV(swest), &
           Vertex(i)%nbrs(seast),Vertex(i)%wgtV(seast), &
           Vertex(i)%nbrs(nwest),Vertex(i)%wgtV(nwest), &
           Vertex(i)%nbrs(neast),Vertex(i)%wgtV(neast)
#endif
    enddo
  98 format(5x,'GRIDVERTEX #',2x,'PART',2x,'DEG',4x,'W',8x,'E',8x, &
		'S',8x,'N',7x,'SW',7x,'SE',7x,'NW',7x,'NE')
  99 format(10x,I3,8x,I1,2x,I1,2x,8(2x,I3,1x,'(',I2,')'))


  end subroutine PrintGridVertex
  subroutine CheckGridNeighbors(Vertex)
  
  implicit none
  type (GridVertex_t), intent(in) :: Vertex(:)

  integer :: i,j,k,nnbrs,inbrs,nvert
  nvert = SIZE(Vertex)

  do i=1,nvert
        nnbrs = SIZE(Vertex(i)%nbrs)
        do j=1,nnbrs
           inbrs = Vertex(i)%nbrs(j)
           if(inbrs > 0) then
              do k=j+1,nnbrs
                 if( inbrs .eq. Vertex(i)%nbrs(k)) &
                    write(iulog,*)'CheckGridNeighbors: ERROR identical neighbors detected  for Vertex ',i
              enddo
           endif
        enddo
  enddo

  end subroutine CheckGridNeighbors


  subroutine initgridedge(GridEdge,GridVertex)
  implicit none

  type (GridEdge_t), intent(inout)       :: GridEdge(:)
  type (GridVertex_t), intent(in),target :: GridVertex(:)

  integer                                :: i,j,k,iptr,wgtV,wgtP
  integer                                :: nelem,nelem_edge,inbr
  logical                                :: Verbose=.FALSE.

  nelem      = SIZE(GridVertex)
  nelem_edge = SIZE(GridEdge)

  GridEdge(:)%reverse=.FALSE.

  iptr=1
  do j=1,nelem
     do i=1,GridVertex(j)%degree
        if(GridVertex(j)%wgtV(i) .gt. 0) then    ! Do this only if has a non-zero weight
           GridEdge(iptr)%tail      => GridVertex(j)
           GridEdge(iptr)%tail_face =  i
#ifdef TESTGRID
           wgtV                     =  GridVertex(j)%wgtV(i)
           GridEdge(iptr)%wgtV      =  wgtV
           wgtP                     =  GridVertex(j)%wgtP(i)
           GridEdge(iptr)%wgtP      =  wgtP
#endif
           inbr                     =  GridVertex(j)%nbrs(i)
           GridEdge(iptr)%head      => GridVertex(inbr)



#ifdef TESTGRID
           ! allocate the indirect addressing arrays
           allocate(GridEdge(iptr)%TailIndex%ixV(wgtV))
           allocate(GridEdge(iptr)%TailIndex%ixP(wgtP))

           allocate(GridEdge(iptr)%TailIndex%iyV(wgtV))
           allocate(GridEdge(iptr)%TailIndex%iyP(wgtP))

           allocate(GridEdge(iptr)%HeadIndex%ixV(wgtV))
           allocate(GridEdge(iptr)%HeadIndex%ixP(wgtP))

           allocate(GridEdge(iptr)%HeadIndex%iyV(wgtV))
           allocate(GridEdge(iptr)%HeadIndex%iyP(wgtP))
#endif
            ! ===========================================
            ! Need this aweful piece of code to determine
            ! which "face" of the neighbor element the
            ! edge links (i.e. the "head_face")
            ! ===========================================
 
            do k=1,GridVertex(inbr)%degree
               if(GridVertex(inbr)%nbrs(k) == GridVertex(j)%number) then
                  GridEdge(iptr)%head_face=k
               endif
            enddo
 
           iptr=iptr+1

        endif
     end do
  end do


  if (Verbose) then

     print *
     write(iulog,*)"element edge tail,head list: (TEST)"
     do i=1,nelem_edge
        write(iulog,*)GridEdge(i)%tail%number,GridEdge(i)%head%number
     end do

     print *
     write(iulog,*)"element edge tail_face, head_face list: (TEST)"
     do i=1,nelem_edge
        write(iulog,*)GridEdge(i)%tail_face,GridEdge(i)%head_face
     end do
  end if

  end subroutine initgridedge

end module GridGraph_mod
