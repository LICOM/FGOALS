#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module edge_mod
  use kinds, only : int_kind, log_kind, real_kind
  use dimensions_mod, only : MaxNeighEdges
  implicit none
  private

  type, public :: rotation_t
     sequence
     integer  :: nbr                                        ! nbr direction: north south east west
     integer  :: reverse                                    ! 0 = do not reverse order
     ! 1 = reverse order
     real (kind=real_kind), dimension(:,:,:), pointer :: R  ! rotation matrix
  end type rotation_t

  type, public :: EdgeDescriptor_t
     sequence
     integer(kind=int_kind)  :: use_rotation
     integer(kind=int_kind)  :: padding
     integer(kind=int_kind)  :: putmapV(MaxNeighEdges)
     integer(kind=int_kind)  :: putmapP(MaxNeighEdges)
     integer(kind=int_kind)  :: getmapV(MaxNeighEdges)
     integer(kind=int_kind)  :: getmapP(MaxNeighEdges)
     logical(kind=log_kind)  :: reverse(MaxNeighEdges)
     type (rotation_t), dimension(:), pointer :: rot ! Identifies list of edges
     !  that must be rotated, and how
  end type EdgeDescriptor_t


  type, public :: EdgeBuffer_t
     sequence
     real (kind=real_kind), dimension(:,:), pointer :: buf
     real (kind=real_kind), dimension(:,:), pointer :: receive
     integer :: nlyr
     integer :: nbuf
  end type EdgeBuffer_t

  type, public :: LongEdgeBuffer_t
     sequence
     integer :: nlyr
     integer :: nbuf
     integer (kind=int_kind), dimension(:,:), pointer :: buf
     integer (kind=int_kind), dimension(:,:), pointer :: receive
  end type LongEdgeBuffer_t

  public :: initEdgeBuffer, initLongEdgeBuffer
  public :: FreeEdgeBuffer, FreeLongEdgeBuffer
  public :: edgeVpack,  edgeDGVpack, LongEdgeVpack


  public :: edgeVunpack,edgeDGVunpack, edgeVunpackVert
  public :: edgeVunpackMIN, LongEdgeVunpackMIN
  public :: edgeVunpackMAX

  public :: edgerotate
  public :: buffermap
  logical, private :: threadsafe=.true.

contains

  ! =========================================
  ! initEdgeBuffer:
  !
  ! create an Real based communication buffer
  ! =========================================
  subroutine initEdgeBuffer(edge,nlyr)
    use dimensions_mod, only : nv, nelemd

    implicit none
    integer,intent(in)                :: nlyr
    type (EdgeBuffer_t),intent(out) :: edge

    ! Local variables

    integer :: nbuf

    nbuf=4*(nv+1)*nelemd
    edge%nlyr=nlyr
    edge%nbuf=nbuf
    allocate(edge%buf(nlyr,nbuf))
    edge%buf(:,:)=0.0D0

    allocate(edge%receive(nlyr,nbuf))
    edge%receive(:,:)=0.0D0

  end subroutine initEdgeBuffer
  ! =========================================
  ! initLongEdgeBuffer:
  !
  ! create an Integer based communication buffer
  ! =========================================
  subroutine initLongEdgeBuffer(edge,nlyr)
    use dimensions_mod, only : nv, nelemd
    implicit none
    integer,intent(in)                :: nlyr
    type (LongEdgeBuffer_t),intent(out) :: edge

    ! Local variables

    integer :: nbuf

    nbuf=4*(nv+1)*nelemd
    edge%nlyr=nlyr
    edge%nbuf=nbuf
    allocate(edge%buf(nlyr,nbuf))
    edge%buf(:,:)=0

    allocate(edge%receive(nlyr,nbuf))
    edge%receive(:,:)=0

  end subroutine initLongEdgeBuffer
  ! =========================================
  ! edgeDGVpack:
  !
  ! Pack edges of v into buf for DG stencil
  ! =========================================
  subroutine edgeDGVpack(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : nv
    type (EdgeBuffer_t)                      :: edge
    integer,              intent(in)   :: vlyr
    real (kind=real_kind),intent(in)   :: v(nv,nv,vlyr)
    integer,              intent(in)   :: kptr
    type (EdgeDescriptor_t)            :: desc

    ! =========================================
    ! This code is just a wrapper call the 
    !   normal edgeVpack
    ! =========================================
    call edgeVpack(edge,v,vlyr,kptr,desc)

  end subroutine edgeDGVpack

  ! ===========================================
  !  FreeEdgeBuffer:
  !
  !  Freed an edge communication buffer
  ! =========================================
  subroutine FreeEdgeBuffer(edge) 
    implicit none
    type (EdgeBuffer_t),intent(inout) :: edge

    edge%nbuf=0
    edge%nlyr=0


    deallocate(edge%buf)
    deallocate(edge%receive)

  end subroutine FreeEdgeBuffer
  ! ===========================================
  !  FreeLongEdgeBuffer:
  !
  !  Freed an edge communication buffer
  ! =========================================
  subroutine FreeLongEdgeBuffer(edge) 
    implicit none
    type (LongEdgeBuffer_t),intent(inout) :: edge

    edge%nbuf=0
    edge%nlyr=0
    deallocate(edge%buf)
    deallocate(edge%receive)

  end subroutine FreeLongEdgeBuffer

  ! =========================================
  ! edgeVpack:
  !
  ! Pack edges of v into buf...
  ! =========================================
  subroutine edgeVpack(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : nv
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest



    type (EdgeBuffer_t)                      :: edge
    integer,              intent(in)   :: vlyr
    real (kind=real_kind),intent(in)   :: v(nv,nv,vlyr)
    integer,              intent(in)   :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc

    ! Local variables
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,k,ir

    integer :: is,ie,in,iw

    if(.not. threadsafe) then
!$OMP BARRIER
       threadsafe=.true.
    end if

    is = desc%putmapV(south)
    ie = desc%putmapV(east)
    in = desc%putmapV(north)
    iw = desc%putmapV(west)

    if(MODULO(nv,2) == 0 .and. UseUnroll) then 
       do k=1,vlyr
          do i=1,nv,2
             edge%buf(kptr+k,is+i)   = v(i  ,1 ,k)
             edge%buf(kptr+k,is+i+1) = v(i+1,1 ,k)
             edge%buf(kptr+k,ie+i)   = v(nv ,i ,k)
             edge%buf(kptr+k,ie+i+1) = v(nv ,i+1 ,k)
             edge%buf(kptr+k,in+i)   = v(i  ,nv,k)
             edge%buf(kptr+k,in+i+1) = v(i+1  ,nv,k)
             edge%buf(kptr+k,iw+i)   = v(1  ,i ,k)
             edge%buf(kptr+k,iw+i+1) = v(1  ,i+1 ,k)

          enddo
       end do
    else
       do k=1,vlyr
          do i=1,nv
             edge%buf(kptr+k,is+i)   = v(i  ,1 ,k)
             edge%buf(kptr+k,ie+i)   = v(nv ,i ,k)
             edge%buf(kptr+k,in+i)   = v(i  ,nv,k)
             edge%buf(kptr+k,iw+i)   = v(1  ,i ,k)
          enddo
       end do
    endif


    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing

    if(desc%reverse(south)) then
       is = desc%putmapV(south)
       do k=1,vlyr
          do i=1,nv
             ir = nv-i+1
             edge%buf(kptr+k,is+ir)=v(i,1,k)
          enddo
       enddo
    endif

    if(desc%reverse(east)) then
       ie = desc%putmapV(east)
       do k=1,vlyr
          do i=1,nv
             ir = nv-i+1
             edge%buf(kptr+k,ie+ir)=v(nv,i,k)
          enddo
       enddo
    endif

    if(desc%reverse(north)) then
       in = desc%putmapV(north)
       do k=1,vlyr
          do i=1,nv
             ir = nv-i+1
             edge%buf(kptr+k,in+ir)=v(i,nv,k)
          enddo
       enddo
    endif

    if(desc%reverse(west)) then
       iw = desc%putmapV(west)
       do k=1,vlyr
          do i=1,nv
             ir = nv-i+1
             edge%buf(kptr+k,iw+ir)=v(1,i,k)
          enddo
       enddo
    endif

    if (desc%putmapV(swest) /= -1) then
       do k=1,vlyr
          edge%buf(kptr+k,desc%putmapV(swest)+1)=v(1  ,1 ,k)
       end do
    end if

    if (desc%putmapV(seast) /= -1) then
       do k=1,vlyr
          edge%buf(kptr+k,desc%putmapV(seast)+1)=v(nv ,1 ,k)
       end do
    end if

    if (desc%putmapV(neast) /= -1) then
       do k=1,vlyr
          edge%buf(kptr+k,desc%putmapV(neast)+1)=v(nv ,nv,k)
       end do
    end if

    if (desc%putmapV(nwest) /= -1) then
       do k=1,vlyr
          edge%buf(kptr+k,desc%putmapV(nwest)+1)=v(1  ,nv,k)
       end do
    end if

  end subroutine edgeVpack
  ! =========================================
  ! LongEdgeVpack:
  !
  ! Pack edges of v into buf...
  ! =========================================
  subroutine LongEdgeVpack(edge,v,vlyr,kptr,desc)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod, only : nv

    type (LongEdgeBuffer_t)            :: edge
    integer,              intent(in)   :: vlyr
    integer (kind=int_kind),intent(in)   :: v(nv,nv,vlyr)
    integer,              intent(in)   :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc

    ! Local variables
    logical, parameter :: UseUnroll = .TRUE.

    integer :: i,k,ir

    integer :: is,ie,in,iw

    if(.not. threadsafe) then
!$OMP BARRIER
       threadsafe=.true.
    end if

    is = desc%putmapV(south)
    ie = desc%putmapV(east)
    in = desc%putmapV(north)
    iw = desc%putmapV(west)

    if(MODULO(nv,2) == 0 .and. UseUnroll) then 
       do k=1,vlyr
          do i=1,nv,2
             edge%buf(kptr+k,is+i)   = v(i  ,1 ,k)
             edge%buf(kptr+k,is+i+1) = v(i+1,1 ,k)
             edge%buf(kptr+k,ie+i)   = v(nv ,i ,k)
             edge%buf(kptr+k,ie+i+1) = v(nv ,i+1 ,k)
             edge%buf(kptr+k,in+i)   = v(i  ,nv,k)
             edge%buf(kptr+k,in+i+1) = v(i+1  ,nv,k)
             edge%buf(kptr+k,iw+i)   = v(1  ,i ,k)
             edge%buf(kptr+k,iw+i+1) = v(1  ,i+1 ,k)

          enddo
       end do
    else
       do k=1,vlyr
          do i=1,nv
             edge%buf(kptr+k,is+i)   = v(i  ,1 ,k)
             edge%buf(kptr+k,ie+i)   = v(nv ,i ,k)
             edge%buf(kptr+k,in+i)   = v(i  ,nv,k)
             edge%buf(kptr+k,iw+i)   = v(1  ,i ,k)
          enddo
       end do

    endif


    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing

    if(desc%reverse(south)) then
       is = desc%putmapV(south)
       do k=1,vlyr
          do i=1,nv
             ir = nv-i+1
             edge%buf(kptr+k,is+ir)=v(i,1,k)
          enddo
       enddo
    endif

    if(desc%reverse(east)) then
       ie = desc%putmapV(east)
       do k=1,vlyr
          do i=1,nv
             ir = nv-i+1
             edge%buf(kptr+k,ie+ir)=v(nv,i,k)
          enddo
       enddo
    endif

    if(desc%reverse(north)) then
       in = desc%putmapV(north)
       do k=1,vlyr
          do i=1,nv
             ir = nv-i+1
             edge%buf(kptr+k,in+ir)=v(i,nv,k)
          enddo
       enddo
    endif

    if(desc%reverse(west)) then
       iw = desc%putmapV(west)
       do k=1,vlyr
          do i=1,nv
             ir = nv-i+1
             edge%buf(kptr+k,iw+ir)=v(1,i,k)
          enddo
       enddo
    endif

    if (desc%putmapV(swest) /= -1) then
       do k=1,vlyr
          edge%buf(kptr+k,desc%putmapV(swest)+1)=v(1  ,1 ,k)
       end do
    end if

    if (desc%putmapV(seast) /= -1) then
       do k=1,vlyr
          edge%buf(kptr+k,desc%putmapV(seast)+1)=v(nv ,1 ,k)
       end do
    end if

    if (desc%putmapV(neast) /= -1) then
       do k=1,vlyr
          edge%buf(kptr+k,desc%putmapV(neast)+1)=v(nv ,nv,k)
       end do
    end if

    if (desc%putmapV(nwest) /= -1) then
       do k=1,vlyr
          edge%buf(kptr+k,desc%putmapV(nwest)+1)=v(1  ,nv,k)
       end do
    end if

  end subroutine LongEdgeVpack

  ! ========================================
  ! edgeVunpack:
  !
  ! Unpack edges from edge buffer into v...
  ! ========================================

  subroutine edgeVunpack(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : nv
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    type (EdgeBuffer_t),         intent(in)  :: edge

    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(inout) :: v(nv,nv,vlyr)
    integer,               intent(in)  :: kptr
    type (EdgeDescriptor_t)            :: desc

    ! Local
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,k
    integer :: is,ie,in,iw

    threadsafe=.false.

    is=desc%getmapV(south)
    ie=desc%getmapV(east)
    in=desc%getmapV(north)
    iw=desc%getmapV(west)

    if(MODULO(nv,2) == 0 .and. UseUnroll) then 
       do k=1,vlyr
          do i=1,nv,2
             v(i  ,1  ,k) = v(i  ,1  ,k)+edge%buf(kptr+k,is+i  )
             v(i+1,1  ,k) = v(i+1,1  ,k)+edge%buf(kptr+k,is+i+1)
             v(nv ,i  ,k) = v(nv ,i  ,k)+edge%buf(kptr+k,ie+i  )
             v(nv ,i+1,k) = v(nv ,i+1,k)+edge%buf(kptr+k,ie+i+1)
             v(i  ,nv ,k) = v(i  ,nv ,k)+edge%buf(kptr+k,in+i  )
             v(i+1,nv ,k) = v(i+1,nv ,k)+edge%buf(kptr+k,in+i+1)
             v(1  ,i  ,k) = v(1  ,i  ,k)+edge%buf(kptr+k,iw+i  )
             v(1  ,i+1,k) = v(1  ,i+1,k)+edge%buf(kptr+k,iw+i+1)
          end do
       end do
    else
       do k=1,vlyr
          do i=1,nv
             v(i  ,1  ,k) = v(i  ,1  ,k)+edge%buf(kptr+k,is+i  )
             v(nv ,i  ,k) = v(nv ,i  ,k)+edge%buf(kptr+k,ie+i  )
             v(i  ,nv ,k) = v(i  ,nv ,k)+edge%buf(kptr+k,in+i  )
             v(1  ,i  ,k) = v(1  ,i  ,k)+edge%buf(kptr+k,iw+i  )
          end do
       end do
    endif

    if(desc%getmapV(swest) /= -1) then 
       do k=1,vlyr
          v(1  ,1 ,k)=v(1 ,1 ,k)+edge%buf(kptr+k,desc%getmapV(swest)+1)
       enddo
    endif

    if(desc%getmapV(seast) /= -1) then 
       do k=1,vlyr
          v(nv ,1 ,k)=v(nv,1 ,k)+edge%buf(kptr+k,desc%getmapV(seast)+1)
       enddo
    endif

    if(desc%getmapV(neast) /= -1) then 
       do k=1,vlyr
          v(nv ,nv,k)=v(nv,nv,k)+edge%buf(kptr+k,desc%getmapV(neast)+1)
       enddo
    endif

    if(desc%getmapV(nwest) /= -1) then 
       do k=1,vlyr
          v(1  ,nv,k)=v(1 ,nv,k)+edge%buf(kptr+k,desc%getmapV(nwest)+1)
       enddo
    endif

  end subroutine edgeVunpack
  subroutine edgeVunpackVert(edge,v,desc)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod, only : nv
    use coordinate_systems_mod, only: cartesian3D_t 

    type (EdgeBuffer_t),   intent(inout)  :: edge
    type (cartesian3D_t), intent(out) :: v(:,:,:)
    type (EdgeDescriptor_t)            :: desc

    ! Local
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,k
    integer :: is,ie,in,iw

    threadsafe=.false.

    is=desc%getmapV(south)
    ie=desc%getmapV(east)
    in=desc%getmapV(north)
    iw=desc%getmapV(west)


    ! N+S
    do i=1,nv/2
       ! North
       v(3,i ,nv)%x = edge%buf(1,in+i) 
       v(3,i ,nv)%y = edge%buf(2,in+i) 
       v(3,i ,nv)%z = edge%buf(3,in+i) 
       ! South
       v(2,i ,1)%x  = edge%buf(1,is+i) 
       v(2,i ,1)%y  = edge%buf(2,is+i) 
       v(2,i ,1)%z  = edge%buf(3,is+i) 
    enddo

    do i=nv/2+1,nv
       ! North
       v(4,i ,nv)%x = edge%buf(1,in+i) 
       v(4,i ,nv)%y = edge%buf(2,in+i) 
       v(4,i ,nv)%z = edge%buf(3,in+i) 
       ! South
       v(1,i ,1)%x  = edge%buf(1,is+i) 
       v(1,i ,1)%y  = edge%buf(2,is+i) 
       v(1,i ,1)%z  = edge%buf(3,is+i)        
    enddo

    do i=1,nv/2
       ! East
       v(3,nv,i)%x = edge%buf(1,ie+i)
       v(3,nv,i)%y = edge%buf(2,ie+i)
       v(3,nv,i)%z = edge%buf(3,ie+i)       
       ! West
       v(4,1,i)%x  = edge%buf(1,iw+i)
       v(4,1,i)%y  = edge%buf(2,iw+i)
       v(4,1,i)%z  = edge%buf(3,iw+i)
    end do

    do i=nv/2+1,nv
       ! East
       v(2,nv,i)%x = edge%buf(1,ie+i)
       v(2,nv,i)%y = edge%buf(2,ie+i)
       v(2,nv,i)%z = edge%buf(3,ie+i)       
       ! West
       v(1,1,i)%x  = edge%buf(1,iw+i)
       v(1,1,i)%y  = edge%buf(2,iw+i)
       v(1,1,i)%z  = edge%buf(3,iw+i)
    end do

    if(desc%getmapV(swest) /= -1) then 
       v(1,1,1)%x=edge%buf(1,desc%getmapV(swest)+1)
       v(1,1,1)%y=edge%buf(2,desc%getmapV(swest)+1)
       v(1,1,1)%z=edge%buf(3,desc%getmapV(swest)+1)
    else
       v(1,1,1)%x=0_real_kind
       v(1,1,1)%y=0_real_kind
       v(1,1,1)%z=0_real_kind
    endif

    if(desc%getmapV(seast) /= -1) then 
       v(2,nv,1)%x=edge%buf(1,desc%getmapV(seast)+1)
       v(2,nv,1)%y=edge%buf(2,desc%getmapV(seast)+1)
       v(2,nv,1)%z=edge%buf(3,desc%getmapV(seast)+1)
    else
       v(2,nv,1)%x=0_real_kind
       v(2,nv,1)%y=0_real_kind
       v(2,nv,1)%z=0_real_kind
    endif

    if(desc%getmapV(neast) /= -1) then 
       v(3,nv,nv)%x=edge%buf(1,desc%getmapV(neast)+1)
       v(3,nv,nv)%y=edge%buf(2,desc%getmapV(neast)+1)
       v(3,nv,nv)%z=edge%buf(3,desc%getmapV(neast)+1)
    else
       v(3,nv,nv)%x=0_real_kind
       v(3,nv,nv)%y=0_real_kind
       v(3,nv,nv)%z=0_real_kind
    endif

    if(desc%getmapV(nwest) /= -1) then 
       v(4,1,nv)%x=edge%buf(1,desc%getmapV(nwest)+1)
       v(4,1,nv)%y=edge%buf(2,desc%getmapV(nwest)+1)
       v(4,1,nv)%z=edge%buf(3,desc%getmapV(nwest)+1)
    else
       v(4,1,nv)%x=0_real_kind
       v(4,1,nv)%y=0_real_kind
       v(4,1,nv)%z=0_real_kind
    endif

    ! Fill the missing vertex info

    do i=2,nv/2
       ! North
       v(4,i ,nv)%x = v(3,i-1 ,nv)%x 
       v(4,i ,nv)%y = v(3,i-1 ,nv)%y
       v(4,i ,nv)%z = v(3,i-1 ,nv)%z
       ! South
       v(1,i ,1)%x  = v(2,i-1 ,1)%x 
       v(1,i ,1)%y  = v(2,i-1 ,1)%y 
       v(1,i ,1)%z  = v(2,i-1 ,1)%z 
    enddo

    do i=nv/2+1,nv-1
       ! North
       v(3,i ,nv)%x = v(4,i+1 ,nv)%x 
       v(3,i ,nv)%y = v(4,i+1 ,nv)%y
       v(3,i ,nv)%z = v(4,i+1 ,nv)%z
       ! South
       v(2,i ,1)%x  = v(1,i+1 ,1)%x 
       v(2,i ,1)%y  = v(1,i+1 ,1)%y
       v(2,i ,1)%z  = v(1,i+1 ,1)%z
    enddo

    do i=2,nv/2
       ! East
       v(2,nv,i)%x = v(3,nv,i-1)%x
       v(2,nv,i)%y = v(3,nv,i-1)%y
       v(2,nv,i)%z = v(3,nv,i-1)%z
       ! West
       v(1,1,i)%x  = v(4,1,i-1)%x
       v(1,1,i)%y  = v(4,1,i-1)%y
       v(1,1,i)%z  = v(4,1,i-1)%z
    end do

    do i=nv/2+1,nv-1
       ! East
       v(3,nv,i)%x = v(2,nv,i+1)%x 
       v(3,nv,i)%y = v(2,nv,i+1)%y
       v(3,nv,i)%z = v(2,nv,i+1)%z
       ! West
       v(4,1,i)%x  = v(1,1,i+1)%x 
       v(4,1,i)%y  = v(1,1,i+1)%y
       v(4,1,i)%z  = v(1,1,i+1)%z
    end do

  end subroutine edgeVunpackVert
  ! ========================================
  ! edgeDGVunpack:
  !
  ! Unpack edges from edge buffer into v...
  ! ========================================

  subroutine edgeDGVunpack(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : nv
    use control_mod, only : north, south, east, west

    type (EdgeBuffer_t),         intent(in)  :: edge
    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(out) :: v(0:nv+1,0:nv+1,vlyr)
    integer,               intent(in)  :: kptr
    type (EdgeDescriptor_t)            :: desc

    ! Local

    integer :: i,k
    integer :: is,ie,in,iw

    threadsafe=.false.

    is=desc%getmapV(south)
    ie=desc%getmapV(east)
    in=desc%getmapV(north)
    iw=desc%getmapV(west)
    do k=1,vlyr
       do i=1,nv
          v(i   ,0   ,k)=edge%buf(kptr+k,is+i)
          v(nv+1,i   ,k)=edge%buf(kptr+k,ie+i)
          v(i   ,nv+1,k)=edge%buf(kptr+k,in+i)
          v(0   ,i   ,k)=edge%buf(kptr+k,iw+i)
       end do
    end do

  end subroutine edgeDGVunpack

  ! ========================================
  ! edgeVunpackMIN/MAX:
  !
  ! Finds the Min/Max edges from edge buffer into v...
  ! ========================================

  subroutine edgeVunpackMAX(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : nv
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t),         intent(in)  :: edge
    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(out) :: v(nv,nv,vlyr)
    integer,               intent(in)  :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc


    ! Local

    integer :: i,k
    integer :: is,ie,in,iw

    threadsafe=.false.

    is=desc%getmapV(south)
    ie=desc%getmapV(east)
    in=desc%getmapV(north)
    iw=desc%getmapV(west)
    do k=1,vlyr
       do i=1,nv
          v(i  ,1  ,k) = MAX(v(i  ,1  ,k),edge%buf(kptr+k,is+i  ))
          v(nv ,i  ,k) = MAX(v(nv ,i  ,k),edge%buf(kptr+k,ie+i  ))
          v(i  ,nv ,k) = MAX(v(i  ,nv ,k),edge%buf(kptr+k,in+i  ))
          v(1  ,i  ,k) = MAX(v(1  ,i  ,k),edge%buf(kptr+k,iw+i  ))
       end do
    end do

    if(desc%getmapV(swest) /= -1) then 
       do k=1,vlyr
          v(1  ,1 ,k)=MAX(v(1 ,1 ,k),edge%buf(kptr+k,desc%getmapV(swest)+1))
       enddo
    endif

    if(desc%getmapV(seast) /= -1) then 
       do k=1,vlyr
          v(nv ,1 ,k)=MAX(v(nv,1 ,k),edge%buf(kptr+k,desc%getmapV(seast)+1))
       enddo
    endif

    if(desc%getmapV(neast) /= -1) then 
       do k=1,vlyr
          v(nv ,nv,k)=MAX(v(nv,nv,k),edge%buf(kptr+k,desc%getmapV(neast)+1))
       enddo
    endif

    if(desc%getmapV(nwest) /= -1) then 
       do k=1,vlyr
          v(1  ,nv,k)=MAX(v(1 ,nv,k),edge%buf(kptr+k,desc%getmapV(nwest)+1))
       enddo
    endif

  end subroutine edgeVunpackMAX
  subroutine edgeVunpackMIN(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : nv
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t),         intent(in)  :: edge
    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(out) :: v(nv,nv,vlyr)
    integer,               intent(in)  :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc


    ! Local

    integer :: i,k
    integer :: is,ie,in,iw

    threadsafe=.false.

    is=desc%getmapV(south)
    ie=desc%getmapV(east)
    in=desc%getmapV(north)
    iw=desc%getmapV(west)
    do k=1,vlyr
       do i=1,nv
          v(i  ,1  ,k) = MIN(v(i  ,1  ,k),edge%buf(kptr+k,is+i  ))
          v(nv ,i  ,k) = MIN(v(nv ,i  ,k),edge%buf(kptr+k,ie+i  ))
          v(i  ,nv ,k) = MIN(v(i  ,nv ,k),edge%buf(kptr+k,in+i  ))
          v(1  ,i  ,k) = MIN(v(1  ,i  ,k),edge%buf(kptr+k,iw+i  ))
       end do
    end do

    if(desc%getmapV(swest) /= -1) then 
       do k=1,vlyr
          v(1  ,1 ,k)=MIN(v(1 ,1 ,k),edge%buf(kptr+k,desc%getmapV(swest)+1))
       enddo
    endif

    if(desc%getmapV(seast) /= -1) then 
       do k=1,vlyr
          v(nv ,1 ,k)=MIN(v(nv,1 ,k),edge%buf(kptr+k,desc%getmapV(seast)+1))
       enddo
    endif

    if(desc%getmapV(neast) /= -1) then 
       do k=1,vlyr
          v(nv ,nv,k)=MIN(v(nv,nv,k),edge%buf(kptr+k,desc%getmapV(neast)+1))
       enddo
    endif

    if(desc%getmapV(nwest) /= -1) then 
       do k=1,vlyr
          v(1  ,nv,k)=MIN(v(1 ,nv,k),edge%buf(kptr+k,desc%getmapV(nwest)+1))
       enddo
    endif

  end subroutine edgeVunpackMIN
  ! ========================================
  ! LongEdgeVunpackMIN:
  !
  ! Finds the Min edges from edge buffer into v...
  ! ========================================

  subroutine LongEdgeVunpackMIN(edge,v,vlyr,kptr,desc)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod, only : nv

    type (LongEdgeBuffer_t),         intent(in)  :: edge
    integer,               intent(in)  :: vlyr
    integer (kind=int_kind), intent(inout) :: v(nv,nv,vlyr)
    integer,               intent(in)  :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc


    ! Local

    integer :: i,k
    integer :: is,ie,in,iw

    threadsafe=.false.

    is=desc%getmapV(south)
    ie=desc%getmapV(east)
    in=desc%getmapV(north)
    iw=desc%getmapV(west)
    do k=1,vlyr
       do i=1,nv
          v(i  ,1  ,k) = MIN(v(i  ,1  ,k),edge%buf(kptr+k,is+i  ))
          v(nv ,i  ,k) = MIN(v(nv ,i  ,k),edge%buf(kptr+k,ie+i  ))
          v(i  ,nv ,k) = MIN(v(i  ,nv ,k),edge%buf(kptr+k,in+i  ))
          v(1  ,i  ,k) = MIN(v(1  ,i  ,k),edge%buf(kptr+k,iw+i  ))
       end do
    end do

    if(desc%getmapV(swest) /= -1) then 
       do k=1,vlyr
          v(1  ,1 ,k)=MIN(v(1 ,1 ,k),edge%buf(kptr+k,desc%getmapV(swest)+1))
       enddo
    endif

    if(desc%getmapV(seast) /= -1) then 
       do k=1,vlyr
          v(nv ,1 ,k)=MIN(v(nv,1 ,k),edge%buf(kptr+k,desc%getmapV(seast)+1))
       enddo
    endif

    if(desc%getmapV(neast) /= -1) then 
       do k=1,vlyr
          v(nv ,nv,k)=MIN(v(nv,nv,k),edge%buf(kptr+k,desc%getmapV(neast)+1))
       enddo
    endif

    if(desc%getmapV(nwest) /= -1) then 
       do k=1,vlyr
          v(1  ,nv,k)=MIN(v(1 ,nv,k),edge%buf(kptr+k,desc%getmapV(nwest)+1))
       enddo
    endif

  end subroutine LongEdgeVunpackMIN

  ! =============================
  ! edgerotate:
  !
  ! Rotate edges in buffer...
  ! =============================

  subroutine edgerotate(edge,vlyr,kptr,desc)
    use dimensions_mod, only : nv
    type (EdgeBuffer_t)           :: edge         ! edge struct
    integer, intent(in)           :: vlyr         ! number of 2d vector fields to rotate
    integer, intent(in)           :: kptr         ! layer pointer into edge buffer
    type (EdgeDescriptor_t)       :: desc

    ! Local variables

    integer :: i,k,k1,k2
    integer :: irot,ia,nbr

    real(kind=real_kind), dimension(:,:,:), pointer :: R
    real(kind=real_kind)  :: tmp1,tmp2

#ifdef _USEASSOCIATED
    if (associated(rot)) then
#else
       if (desc%use_rotation == 1) then
#endif

          do irot=1,SIZE(desc%rot)

             nbr  =  desc%rot(irot)%nbr
             R    => desc%rot(irot)%R

             ia=desc%putmapV(nbr)

             ! ========================================
             ! If nbr direction is (1-4) => is an edge
             ! ========================================

             if (nbr <= 4) then

                ! ========================================================
                !  Is an edge. Rotate it in place
                ! ========================================================

                do i=1,nv
                   do k=1,vlyr,2
                      k1 = kptr + k
                      k2 = k1 + 1
                      tmp1=R(1,1,i)*edge%buf(k1,ia+i) + R(1,2,i)*edge%buf(k2,ia+i)
                      tmp2=R(2,1,i)*edge%buf(k1,ia+i) + R(2,2,i)*edge%buf(k2,ia+i)
                      edge%buf(k1,ia+i)=tmp1
                      edge%buf(k2,ia+i)=tmp2
                   end do
                end do

             else

                ! ===================================================
                ! Is an element corner point, but not a cube corner
                ! point, just rotate it in place.
                ! ===================================================

                if (ia /= -1) then
                   do k=1,vlyr,2
                      k1 = kptr + k
                      k2 = k1+1
                      tmp1=R(1,1,1)*edge%buf(k1,ia+1) + R(1,2,1)*edge%buf(k2,ia+1)
                      tmp2=R(2,1,1)*edge%buf(k1,ia+1) + R(2,2,1)*edge%buf(k2,ia+1)
                      edge%buf(k1,ia+1)=tmp1
                      edge%buf(k2,ia+1)=tmp2
                   end do
                end if

             end if

          end do

       endif

     end subroutine edgerotate

     ! =============================================
     ! buffermap:
     !
     ! buffermap translates element number, inum and
     ! element edge/corner, facet, into an edge buffer 
     ! memory location, loc.
     ! =============================================

     function buffermap(inum,facet) result(loc)
       use dimensions_mod, only : nv
       integer, intent(in) :: inum   
       integer, intent(in) :: facet
       integer :: loc

       if (facet>4) then
          if (inum == -1) then
             loc = inum
          else
             loc=(inum-1)*(4*nv+4)+4*nv+(facet-5)
          end if
       else
          loc=(inum-1)*(4*nv+4)+nv*(facet-1)
       end if

     end function buffermap

   end module edge_mod









