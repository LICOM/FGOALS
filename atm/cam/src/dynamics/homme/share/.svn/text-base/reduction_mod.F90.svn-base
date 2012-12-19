#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module reduction_mod
  use kinds, only : real_kind
  implicit none
  private

  type, public :: ReductionBuffer_r_1d_t
     sequence
     real (kind=real_kind), dimension(:), pointer :: buf
     integer :: len
     integer :: ctr
  end type ReductionBuffer_r_1d_t

  type, public :: ReductionBuffer_ordered_1d_t
     sequence
     real (kind=real_kind), dimension(:,:),pointer :: buf
     integer :: len
     integer :: ctr
  end type ReductionBuffer_ordered_1d_t

  public :: ParallelMin,ParallelSum,ParallelMax

  !type (ReductionBuffer_ordered_1d_t), public :: red_sum
  type (ReductionBuffer_r_1d_t),       public :: red_sum
  type (ReductionBuffer_r_1d_t),       public :: red_max,red_min
  type (ReductionBuffer_r_1d_t),       public :: red_flops,red_timer

  !JMD new addition
#ifndef Darwin
  SAVE red_sum,red_max,red_min,red_flops,red_timer
#endif
  interface ParallelMin
     module procedure ParallelMin1d
     module procedure ParallelMin0d
  end interface
  interface ParallelMax
     module procedure ParallelMax1d
     module procedure ParallelMax0d
  end interface
  interface ParallelSum
     module procedure ParallelSum1d
     module procedure ParallelSum0dv
     module procedure ParallelSum0d
  end interface
  interface psum_mt
     module procedure psum_mt_r_1d
     module procedure psum_mt_ordered_1d
!     module procedure psum_mt_r_1d_time
!     module procedure psum_mt_ordered_1d_time
  end interface

  interface pmax_mt
     module procedure pmax_mt_r_1d
  end interface

  interface pmin_mt
     module procedure pmin_mt_r_1d
  end interface

  interface InitReductionBuffer
     module procedure InitReductionBuffer_r_1d
     module procedure InitReductionBuffer_ordered_1d
  end interface

  public :: initReductionBuffer
  public :: psum_mt, pmax_mt, pmin_mt
  public :: ElementSum_1d

contains

  function ParallelMin1d(data,hybrid) result(pmin)
    use hybrid_mod, only : hybrid_t
    implicit none
    real(kind=real_kind), intent(in)    :: data(:)
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=real_kind)                :: pmin

    real(kind=real_kind)                :: tmp(1)


    tmp(1) = MINVAL(data)
    call pmin_mt(red_min,tmp,1,hybrid)
    pmin = red_min%buf(1)

  end function ParallelMin1d

  function ParallelMin0d(data,hybrid) result(pmin)
    use hybrid_mod, only : hybrid_t
    implicit none
    real(kind=real_kind), intent(in)    :: data
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=real_kind)                :: pmin
    real(kind=real_kind)                :: tmp(1)
    tmp(1) = data
    call pmin_mt(red_min,tmp,1,hybrid)
    pmin = red_min%buf(1)

  end function ParallelMin0d
  !==================================================
  function ParallelMax1d(data,hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    implicit none
    real(kind=real_kind), intent(in)    :: data(:)
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=real_kind)                :: pmax

    real(kind=real_kind)                :: tmp(1)


    tmp(1) = MAXVAL(data)
    call pmax_mt(red_max,tmp,1,hybrid)
    pmax = red_max%buf(1)

  end function ParallelMax1d
  function ParallelMax0d(data,hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    implicit none
    real(kind=real_kind), intent(in)    :: data
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=real_kind)                :: pmax
    real(kind=real_kind)                :: tmp(1)

    tmp(1)=data

    call pmax_mt(red_max,tmp,1,hybrid)
    pmax = red_max%buf(1)

  end function ParallelMax0d
  !==================================================
  ! function ParallelSum(data,hybrid) result(psum)
  !   implicit none
  !   real(kind=real_kind), intent(in)    :: data(:)
  !   type (hybrid_t),      intent(in)    :: hybrid
  !   real(kind=real_kind)                :: psum

  !   real(kind=real_kind)                :: tmp(1)


  !  tmp(1) = SUM(data)
  !  call psum_mt_ordered_1d(red_sum,tmp,1,hybrid)
  !  psum = red_sum%buf(1,1)

  ! end function ParallelSum

  !==================================================
  function ParallelSum0dv(data,len,hybrid) result(psum)
    use hybrid_mod, only : hybrid_t
    implicit none
    real(kind=real_kind), intent(in)    :: data(len)
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=real_kind)                :: psum(len),tmp(len)
    integer                             :: len

    ! need to check that size of red_sum%buf >= len
    tmp = data
    call psum_mt(red_sum,tmp,len,hybrid)
    psum(1:len) = red_sum%buf(1:len)

  end function ParallelSum0dv
  function ParallelSum1d(data,hybrid) result(psum)
    use hybrid_mod, only : hybrid_t
    implicit none
    real(kind=real_kind), intent(in)    :: data(:)
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=real_kind)                :: psum

    real(kind=real_kind)                :: tmp(1)


    tmp(1) = SUM(data)
    call psum_mt(red_sum,tmp,1,hybrid)
    psum = red_sum%buf(1)

  end function ParallelSum1d
  function ParallelSum0d(data,hybrid) result(psum)
    use hybrid_mod, only : hybrid_t
    implicit none
    real(kind=real_kind), intent(in)    :: data
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=real_kind)                :: psum
    real(kind=real_kind)                :: tmp(1)

    tmp(1)=data

    call psum_mt(red_sum,tmp,1,hybrid)
    psum = red_sum%buf(1)

  end function ParallelSum0d
  ! ================================
  subroutine initReductionBuffer_r_1d(red,len)
    integer, intent(in)           :: len
    type (ReductionBuffer_r_1d_t),intent(out) :: red
    red%len  = len
    allocate(red%buf(len))
    red%buf  = 0.0D0
    red%ctr  = 0
  end subroutine initReductionBuffer_r_1d
  !****************************************************************
  subroutine initReductionBuffer_ordered_1d(red,len,nthread)
    integer, intent(in)           :: len
    integer, intent(in)           :: nthread
    type (ReductionBuffer_ordered_1d_t),intent(out) :: red
    red%len  = len
    allocate(red%buf(len,nthread+1))
    red%buf  = 0.0D0
    red%ctr  = 0
  end subroutine initReductionBuffer_ordered_1d

  ! =======================================
  ! psum_mt:
  !
  ! thread safe, parallel reduce sum
  ! of a one dimensional reduction vector
  ! =======================================
#ifdef TIMERS
  subroutine psum_mt_r_1d_time(red,redp,len,hybrid)
    use hybrid_mod, only : hybrid_t
#ifdef _MPI
    use parallel_mod, only: mpi_sum, mpireal_t, mpi_min, mpi_max
#endif
    type (ReductionBuffer_r_1d_t),intent(inout)  :: red     ! shared memory reduction buffer struct
    real (kind=real_kind), intent(inout)            :: redp(:) ! thread private vector of partial sum 
    integer,               intent(in)            :: len     ! buffer length
    type (hybrid_t),       intent(in)            :: hybrid  ! parallel handle

#ifdef _MPI
    integer                 :: i
    integer ierr
#endif

    ! Local variables

    integer k

    ! insert length safety check here
    ! compare SIZE(redp) with red%len


    !$OMP BARRIER
    !$OMP CRITICAL (CRITSUM_TIME)
    if (red%ctr == 0) then
       red%buf(1:len)=redp(1:len)
    else
       do k=1,len
          red%buf(k)=red%buf(k)+redp(k)
       end do
    end if
    red%ctr=red%ctr+1
    if (red%ctr == hybrid%NThreads) red%ctr=0
    !$OMP END CRITICAL (CRITSUM_TIME)

    !$OMP BARRIER

#ifdef _MPI
    if (hybrid%ithr==0) then

       call MPI_Allreduce(red%buf(1),redp,len,MPIreal_t, &
            MPI_SUM,hybrid%par%comm,ierr)

       do i=1,len
          red%buf(i)=redp(i)
       enddo
    end if
    !$OMP BARRIER
#endif   

  end subroutine psum_mt_r_1d_time

  !
  subroutine psum_mt_ordered_1d_time(red,redp,len,hybrid)
    use hybrid_mod, only : hybrid_t
#ifdef _MPI
    use parallel_mod, only: mpi_sum, mpireal_t, mpi_min, mpi_max
#endif
    implicit none
    type (ReductionBuffer_ordered_1d_t),intent(inout)  :: red     ! shared memory reduction buffer struct
    real (kind=real_kind), intent(inout)            :: redp(:) ! thread private vector of partial sum 
    integer,               intent(in)            :: len     ! buffer length
    type (hybrid_t),       intent(in)            :: hybrid  ! parallel handle
    type (timer_t),        intent(inout)         :: timer
    !    real*8,                intent(inout)         :: timer(:)
    real*8  :: et,st

    ! Local variables

    integer i
    integer k,slot,nthreads
#ifdef _MPI
    integer ierr
#endif

    ! insert length safety check here
    ! compare SIZE(redp) with red%len

    !$OMP BARRIER
    TIMER_DETAIL_START(timer,1,st)
    !$OMP CRITICAL
    slot = hybrid%ithr+2
    nthreads = hybrid%NThreads
    do k=1,len
       red%buf(k,slot)=redp(k)
    enddo
    red%ctr=red%ctr+1
    if (red%ctr == nthreads) then 
       red%ctr=0
       do k=1,len
          red%buf(k,1) = 0.D0
       enddo
       do i=1,nthreads
          do k=1,len
             red%buf(k,1) = red%buf(k,1) + red%buf(k,i+1)
          enddo
       enddo
    endif
    !$OMP END CRITICAL

    !$OMP BARRIER
    if(TIMER_DETAIL(1,timer)) then
       TIMER_START(et)
       TIMER_UPDATE(timer,T_REDUCE_SMP,et,st)
    endif
#ifdef _MPI
    if (hybrid%ithr==0) then

       call MPI_Allreduce(red%buf(1,1),redp,len,MPIreal_t, &
            MPI_SUM,hybrid%par%comm,ierr)

       red%buf(1:len,1)=redp(1:len)
    end if
    !$OMP BARRIER
#endif   

  end subroutine psum_mt_ordered_1d_time
#endif
  subroutine psum_mt_r_1d(red,redp,len,hybrid)
    use hybrid_mod, only : hybrid_t
#ifdef _MPI
    use parallel_mod, only: mpi_sum, mpireal_t
#endif
    implicit none
    type (ReductionBuffer_r_1d_t),intent(inout)  :: red     ! shared memory reduction buffer struct
    real (kind=real_kind), intent(inout)            :: redp(:) ! thread private vector of partial sum 
    integer,               intent(in)            :: len     ! buffer length
    type (hybrid_t),       intent(in)            :: hybrid  ! parallel handle

    ! Local variables

    integer k
#ifdef _MPI
    integer ierr
#endif

    ! insert length safety check here
    ! compare SIZE(redp) with red%len


    !$OMP BARRIER
    !$OMP CRITICAL (CRITSUM)
    if (red%ctr == 0) then
       red%buf(1:len)=redp(1:len)
    else
       red%buf(1:len)=red%buf(1:len)+redp(1:len)
    end if
    red%ctr=red%ctr+1
    if (red%ctr == hybrid%NThreads) red%ctr=0
    !$OMP END CRITICAL (CRITSUM)

#ifdef _MPI
    !$OMP BARRIER
    if (hybrid%ithr==0) then
       call MPI_Allreduce(red%buf(1),redp,len,MPIreal_t, &
            MPI_SUM,hybrid%par%comm,ierr)

       red%buf(1:len)=redp(1:len)
    end if
#endif   
    !$OMP BARRIER
  end subroutine psum_mt_r_1d

  subroutine psum_mt_ordered_1d(red,redp,len,hybrid)
    use hybrid_mod, only : hybrid_t
#ifdef _MPI
    use parallel_mod, only: mpi_sum, mpireal_t
#endif

    implicit none
    type (ReductionBuffer_ordered_1d_t),intent(inout)  :: red     ! shared memory reduction buffer struct
    real (kind=real_kind), intent(inout)            :: redp(:) ! thread private vector of partial sum 
    integer,               intent(in)            :: len     ! buffer length
    type (hybrid_t),       intent(in)            :: hybrid  ! parallel handle

    ! Local variables

    integer i
    integer k,slot,nthreads
#ifdef _MPI
    integer ierr
#endif

    ! insert length safety check here
    ! compare SIZE(redp) with red%len

    !JMD!$OMP ORDERED
    !JMD!$OMP END ORDERED

    !$OMP BARRIER
    !$OMP CRITICAL
    slot = hybrid%ithr+2
    nthreads = hybrid%NThreads
    do k=1,len
       red%buf(k,slot)=redp(k)
    enddo
    red%ctr=red%ctr+1
    if (red%ctr == nthreads) then 
       red%ctr=0
       do k=1,len
          red%buf(k,1) = 0.D0
       enddo
       do i=1,nthreads
          do k=1,len
             red%buf(k,1) = red%buf(k,1) + red%buf(k,i+1)
          enddo
       enddo
    endif
    !$OMP END CRITICAL

#ifdef _MPI
    !$OMP BARRIER
    if (hybrid%ithr==0) then
       call MPI_Allreduce(red%buf(1,1),redp,len,MPIreal_t, &
            MPI_SUM,hybrid%par%comm,ierr)
       red%buf(1:len,1)=redp(1:len)
    end if
#endif   
    !$OMP BARRIER
  end subroutine psum_mt_ordered_1d
  ! =======================================
  ! pmax_mt:
  !
  ! thread safe, parallel reduce maximum
  ! of a one dimensional reduction vector
  ! =======================================

  subroutine pmax_mt_r_1d(red,redp,len,hybrid)
    use hybrid_mod, only : hybrid_t
#ifdef _MPI
    use parallel_mod, only: mpi_min, mpi_max, mpireal_t
#endif

    type (ReductionBuffer_r_1d_t)     :: red     ! shared memory reduction buffer struct
    real (kind=real_kind), intent(inout) :: redp(:) ! thread private vector of partial sum
    integer,               intent(in) :: len     ! buffer length
    type (hybrid_t),       intent(in) :: hybrid  ! parallel handle

    ! Local variables
#ifdef _MPI
    integer ierr
#endif

    integer  :: k
    ! insert length safety check here
    ! compare SIZE(redp) with red%len

    !$OMP BARRIER
    !$OMP CRITICAL (CRITMAX)
    if (red%ctr == 0) red%buf(1:len)= -9.11e30
    if (red%ctr < hybrid%NThreads) then
       do k=1,len
          red%buf(k)=MAX(red%buf(k),redp(k))
       enddo
       red%ctr=red%ctr+1
    end if
    if (red%ctr == hybrid%NThreads) red%ctr=0
    !$OMP END CRITICAL (CRITMAX)
#ifdef _MPI
    !$OMP BARRIER
    if (hybrid%ithr==0) then

       call MPI_Allreduce(red%buf(1),redp,len,MPIreal_t, &
            MPI_MAX,hybrid%par%comm,ierr)

       red%buf(1:len)=redp(1:len)
    end if
#endif
    !$OMP BARRIER


  end subroutine pmax_mt_r_1d

  ! =======================================
  ! pmin_mt:
  !
  ! thread safe, parallel reduce maximum
  ! of a one dimensional reduction vector
  ! =======================================

  subroutine pmin_mt_r_1d(red,redp,len,hybrid)
    use kinds, only : int_kind
    use hybrid_mod, only : hybrid_t
#ifdef _MPI
    use parallel_mod, only: mpi_min, mpireal_t
#endif

    type (ReductionBuffer_r_1d_t)     :: red     ! shared memory reduction buffer struct
    real (kind=real_kind), intent(inout) :: redp(:) ! thread private vector of partial sum
    integer,               intent(in) :: len     ! buffer length
    type (hybrid_t),       intent(in) :: hybrid  ! parallel handle

    ! Local variables

#ifdef _MPI
    integer ierr
#endif
    integer (kind=int_kind) :: k

    ! insert length safety check here
    ! compare SIZE(redp) with red%len

    !$OMP BARRIER
    !$OMP CRITICAL (CRITMAX)
    if (red%ctr == 0) red%buf(1:len)= 9.11e30
    if (red%ctr < hybrid%NThreads) then
       do k=1,len
          red%buf(k)=MIN(red%buf(k),redp(k))
       enddo
       red%ctr=red%ctr+1
    end if
    if (red%ctr == hybrid%NThreads) red%ctr=0
    !$OMP END CRITICAL (CRITMAX)
#ifdef _MPI
    !$OMP BARRIER
    if (hybrid%ithr==0) then

       call MPI_Allreduce(red%buf(1),redp,len,MPIreal_t, &
            MPI_MIN,hybrid%par%comm,ierr)

       red%buf(1:len)=redp(1:len)
    end if
#endif
    !$OMP BARRIER


  end subroutine pmin_mt_r_1d

  subroutine ElementSum_1d(res,variable,type,hybrid)
    use hybrid_mod, only : hybrid_t
    use dimensions_mod, only : nelem
#ifdef _MPI
  use parallel_mod, only : ORDERED, mpireal_t, mpi_min, mpi_max, mpi_sum, mpi_success
#else
  use parallel_mod, only : ORDERED
#endif
    implicit none

    ! ==========================
    !     Arguments
    ! ==========================
    real(kind=real_kind),intent(out) :: res
    real(kind=real_kind),intent(in)  :: variable(:)
    integer,intent(in)               :: type
    type (hybrid_t), intent(in)      :: hybrid 

    ! ==========================
    !       Local Variables
    ! ==========================

    !
    ! Note this is a real kludge here since it may be used for
    !  arrays of size other then nelem
    !

    integer                          :: i
#if 0
    real(kind=real_kind),allocatable :: Global(:)
    real(kind=real_kind),allocatable :: buffer(:)
#endif

#ifdef _MPI
    integer                           :: errorcode,errorlen
    character*(80) errorstring

    real(kind=real_kind)             :: local_sum
    integer                          :: ierr
#endif

#ifdef _MPI
    if(hybrid%ithr == 0) then 
#if 0
       if(type == ORDERED) then
          allocate(buffer(nelem))
          call MPI_Gatherv(variable,nelemd,MPIreal_t,buffer, &
               recvcount,displs,MPIreal_t,hybrid%par%root, &
               hybrid%par%comm,ierr)
          if(ierr .ne. MPI_SUCCESS) then 
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'ElementSum_1d: Effor after call to MPI_Gatherv: ',errorstring
          endif

          if(hybrid%par%masterproc) then
             allocate(Global(nelem))
             do ip=1,hybrid%par%nprocs
                nelemr = recvcount(ip)
                disp   = displs(ip)
                do ie=1,nelemr
                   ig = Schedule(ip)%Local2Global(ie)
                   Global(ig) = buffer(disp+ie)
                enddo
             enddo
             ! ===========================
             !  Perform the ordererd sum
             ! ===========================
             res = 0.0d0
             do i=1,nelem
                res = res + Global(i)
             enddo
             deallocate(Global)
          endif
          ! =============================================
          !  Broadcast the results back everybody
          ! =============================================
          call MPI_Bcast(res,1,MPIreal_t,hybrid%par%root, &
               hybrid%par%comm,ierr)
          if(ierr .ne. MPI_SUCCESS) then 
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'ElementSum_1d: Effor after call to MPI_Bcast: ',errorstring
          endif

          deallocate(buffer)
       else
#endif
          local_sum=SUM(variable)
          call MPI_Barrier(hybrid%par%comm,ierr)

          call MPI_Allreduce(local_sum,res,1,MPIreal_t, &
               MPI_SUM,hybrid%par%comm,ierr)
          if(ierr .ne. MPI_SUCCESS) then 
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'ElementSum_1d: Effor after call to MPI_Allreduce: ',errorstring
          endif
#if 0
       endif
#endif
    endif
#else
    if(hybrid%ithr == 0) then 
       if(type == ORDERED) then
          ! ===========================
          !  Perform the ordererd sum
          ! ===========================
          res = 0.0d0
          do i=1,nelem
             res = res + variable(i)
          enddo
       else
          res=SUM(variable)
       endif
    endif
#endif

  end subroutine ElementSum_1d

end module reduction_mod
