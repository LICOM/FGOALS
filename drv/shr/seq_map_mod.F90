module seq_map_mod

!---------------------------------------------------------------------
!
! Purpose:
!
! General self normalizing mapping routine with optional fraction
!       
! Author: T. Craig, Jan-28-2011
!
!---------------------------------------------------------------------

  use shr_kind_mod      ,only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_kind_mod      ,only: CL => SHR_KIND_CL, CX => SHR_KIND_CX
  use shr_sys_mod
  use shr_const_mod
  use mct_mod

  use seq_comm_mct, only : logunit, loglevel

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: seq_map_avNorm

  interface seq_map_avNorm ; module procedure &
    seq_map_avNormArr, &
    seq_map_avNormAvF
  end interface

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

#ifdef CPP_VECTOR
  logical :: usevector = .true.
#else
  logical :: usevector = .false.
#endif
  
#ifdef SYSUNICOS
  logical :: usealltoall = .true.
#else
  logical :: usealltoall = .false.
#endif

!=======================================================================
contains
!=======================================================================

!=======================================================================

  subroutine seq_map_avNormAvF(av_i, av_o, sMatP, avf_i, avfifld, avf_o, avfofld, rList)

    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(mct_aVect) , intent(in)    :: av_i  ! input 
    type(mct_aVect) , intent(inout) :: av_o  ! output
    type(mct_sMatp) , intent(inout) :: sMatp ! sMat
    type(mct_aVect) , intent(in)    :: avf_i  ! extra src "weight"
    character(len=*), intent(in)    :: avfifld ! field name in avf_i
    type(mct_aVect) , intent(in),optional :: avf_o   ! extra dst "weight"
    character(len=*), intent(in),optional :: avfofld ! field name in avf_o
    character(len=*), intent(in),optional :: rList ! fields list
    !
    integer(in) :: lsize_i, lsize_f, lsize_o, kf, j
    real(r8),allocatable :: frac_i(:),frac_o(:)
    character(*),parameter :: subName = '(seq_map_avNormArr) '
    !-----------------------------------------------------

    lsize_i = mct_aVect_lsize(av_i)
    lsize_f = mct_aVect_lsize(avf_i)

    if (lsize_i /= lsize_f) then
       write(logunit,*) subname,' ERROR: lsize_i ne lsize_f ',lsize_i,lsize_f
       call shr_sys_abort(subname//' ERROR size_i ne lsize_f')
    endif

    !--- extract frac_i field from avf_i to pass to seq_map_avNormArr ---
    allocate(frac_i(lsize_i))
    do j = 1,lsize_i
       kf = mct_aVect_indexRA(avf_i,trim(avfifld))
       frac_i(j) = avf_i%rAttr(kf,j)
    enddo

    if (present(avf_o)) then

       if (.not.present(avfofld)) then
          write(logunit,*) subname,' ERROR: avfofld with avf_o '
          call shr_sys_abort(subname//' ERROR avfofld with avf_o')
       endif

       lsize_o = mct_aVect_lsize(av_o)
       lsize_f = mct_aVect_lsize(avf_o)

       if (lsize_o /= lsize_f) then
          write(logunit,*) subname,' ERROR: lsize_o ne lsize_f ',lsize_o,lsize_f
          call shr_sys_abort(subname//' ERROR size_o ne lsize_f')
       endif

       allocate(frac_o(lsize_o))
       do j = 1,lsize_o
          kf = mct_aVect_indexRA(avf_o,trim(avfofld))
          frac_o(j) = avf_o%rAttr(kf,j)
       enddo

    endif

    if (present(rList)) then
       if (present(avf_o)) then
          call seq_map_avNormArr(av_i, av_o, sMatP, frac_i, frac_o, rList=rList)
       else
          call seq_map_avNormArr(av_i, av_o, sMatP, frac_i, rList=rList)
       endif
    else
       if (present(avf_o)) then
          call seq_map_avNormArr(av_i, av_o, sMatP, frac_i, frac_o)
       else
          call seq_map_avNormArr(av_i, av_o, sMatP, frac_i)
       endif
    endif

    if (present(avf_o)) then
       deallocate(frac_o)
    end if
    deallocate(frac_i)

  end subroutine seq_map_avNormAvF

!=======================================================================

  subroutine seq_map_avNormArr(av_i, av_o, sMatP, norm_i, norm_o, rList, donorm)

    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(mct_aVect) , intent(in)    :: av_i  ! input 
    type(mct_aVect) , intent(inout) :: av_o  ! output
    type(mct_sMatp) , intent(inout) :: sMatp ! sMat
    real(r8)        , intent(in), optional :: norm_i(:)  ! source "weight"
    real(r8)        , intent(in), optional :: norm_o(:)  ! dst "weight"
    character(len=*), intent(in), optional :: rList ! fields list
    logical         , intent(in), optional :: donorm  ! normalize at end
    !
    ! Local variables
    !
    type(mct_aVect)        :: avp_i , avp_o
    integer(IN)            :: i,j,ier,kf
    integer(IN)            :: lsize_i,lsize_o
    real(r8)               :: norm
    character(CX)          :: lrList
    logical                :: ldonorm
    character(*),parameter :: subName = '(seq_map_avNormArr) '
    character(len=*),parameter :: ffld = 'norm8wt'  ! want something unique
    !-----------------------------------------------------

    lsize_i = mct_aVect_lsize(av_i)
    lsize_o = mct_aVect_lsize(av_o)

    ldonorm = .true.
    if (present(donorm)) then
       ldonorm = donorm
    endif

    if (present(norm_i) .and..not.ldonorm) then
       write(logunit,*) subname,' ERROR: norm_i and donorm = false'
       call shr_sys_abort(subname//' ERROR norm_i and donorm = false')
    endif

    if (present(norm_o) .and. .not.present(norm_i)) then
       write(logunit,*) subname,' ERROR: norm_i must be with norm_o'
       call shr_sys_abort(subname//' ERROR norm_i must be with norm_o')
    endif

    if (present(norm_i)) then
       if (size(norm_i) /= lsize_i) then
          write(logunit,*) subname,' ERROR: size(norm_i) ne lsize_i ',size(norm_i),lsize_i
          call shr_sys_abort(subname//' ERROR size(norm_i) ne lsize_i')
       endif
    endif

    if (present(norm_o)) then
       if (size(norm_o) /= lsize_o) then
          write(logunit,*) subname,' ERROR: size(norm_o) ne lsize_o ',size(norm_o),lsize_o
          call shr_sys_abort(subname//' ERROR size(norm_o) ne lsize_o')
       endif
    endif

    !--- create temporary avs for mapping ---

    if (present(rList)) then
       call mct_aVect_init(avp_i, rList=trim( rList)//':'//ffld, lsize=lsize_i)
       call mct_aVect_init(avp_o, rList=trim( rList)//':'//ffld, lsize=lsize_o)
    else
       lrList = trim(mct_aVect_exportRList2c(av_i))
       call mct_aVect_init(avp_i, rList=trim(lrList)//':'//ffld, lsize=lsize_i)
       lrList = trim(mct_aVect_exportRList2c(av_o))
       call mct_aVect_init(avp_o, rList=trim(lrList)//':'//ffld, lsize=lsize_o)
    endif

    !--- copy av_i to avp_i and set ffld value to 1.0
    !--- then multiply all fields by norm_i if norm_i exists 
    !--- this will do the right thing for the norm_i normalization 

    call mct_aVect_copy(aVin=av_i, aVout=avp_i, VECTOR=usevector)
    kf = mct_aVect_indexRA(avp_i,ffld)
    do j = 1,lsize_i
       avp_i%rAttr(kf,j) = 1.0_r8
    enddo
    if (present(norm_i)) then
       do j = 1,lsize_i
          avp_i%rAttr(:,j) = avp_i%rAttr(:,j)*norm_i(j)
       enddo
    endif

    !--- map ---

    call mct_sMat_avMult(avp_i, sMatp, avp_o, VECTOR=usevector)

    !--- renormalize avp_o by mapped norm_i, which could be all 1s ---

    if (ldonorm) then
    do j = 1,lsize_o
       if (present(norm_o)) then
          norm = norm_o(j)
       else
          kf = mct_aVect_indexRA(avp_o,ffld)
          norm = avp_o%rAttr(kf,j)
       endif
       if (norm /= 0.0_r8) then
          norm = 1.0_r8/norm
       endif
       avp_o%rAttr(:,j) = avp_o%rAttr(:,j)*norm
    enddo
    endif

    !--- copy back into av_o and we are done ---

    call mct_aVect_copy(aVin=avp_o, aVout=av_o, VECTOR=usevector)

    call mct_aVect_clean(avp_i)
    call mct_aVect_clean(avp_o)

 end subroutine seq_map_avNormArr

!===============================================================================

end module seq_map_mod
