module mct2esmf_mod
#ifdef USE_ESMF_LIB

use ESMF_Mod
use esmfshr_util_mod,               only : esmfshr_util_ArrayGetSize, esmfshr_util_ArrayCopy
use mct_mod
use esmfshr_util_mod
use shr_string_mod, only: shr_string_listGetNum
use shr_log_mod, only: shr_log_unit, shr_log_level
use shr_sys_mod, only: shr_sys_flush, shr_sys_abort

implicit none

! 
! Author: Fei Liu
!
! This module implements methods that work with ESMF_Array and MCT attribute
! vectors (attrVect) and global segmentation map (gsMap)
! These methods converts attrVect and gsMap to Array and distgrid.
! Another layer of utils take these 
! array related methods and put them into interface blocks. user using
! the higher level utils module see/use names such as esmf2mct_init, etc

public mct2esmf_copy
public mct2esmf_init

private

interface mct2esmf_init
    module procedure mct2esmf_init_Distgrid_GSmap
    module procedure mct2esmf_init_Distgrid_Gindex
    module procedure mct2esmf_init_Array_AvGSmap
    module procedure mct2esmf_init_Array_ArrayGSmap
    module procedure mct2esmf_init_Arrays_AvsGSmap
    module procedure mct2esmf_init_Array_GindexList
    module procedure mct2esmf_init_Array_DistgridList
    module procedure mct2esmf_init_Array_gsMapNF
    module procedure mct2esmf_init_Array_gsMaprList
end interface

interface mct2esmf_copy
    module procedure mct2esmf_copy_Array_Av
    module procedure mct2esmf_copy_alist
end interface

contains
!--------------------------------------------------------------------------

! create a ESMF_DistGrid from a gsMap
function mct2esmf_init_Distgrid_GSmap(gsMap, rc)

    implicit none

    ! inout parameters
    type(mct_gsMap), intent(in)             :: gsMap
    integer, intent(out), optional          :: rc

    ! return
    type(ESMF_DistGrid)                     :: mct2esmf_init_Distgrid_GSmap

    ! local variables
    integer                                 :: localrc, listsize, lpet
    integer                                 :: vm_comm_id, mct_comm_id
    type(ESMF_VM)                           :: vm
    integer, pointer                        :: gindex(:)
    character(len=*),parameter :: subname = 'mct2esmf_init_Distgrid_GSmap'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    call ESMF_VMGet(vm, localPet=lpet, mpiCommunicator=vm_comm_id, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    call MPI_Comm_Dup(vm_comm_id, mct_comm_id, localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    listsize = mct_gsMap_lsize(gsMap, comm=mct_comm_id)
    allocate(gindex(listsize), stat=localrc)
    if(localrc /= 0) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)
    call mct_gsMap_OrderedPoints(gsMap, PEno=lpet, Points=gindex)

    mct2esmf_init_Distgrid_GSmap = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    deallocate(gindex)

end function mct2esmf_init_Distgrid_GSmap

!--------------------------------------------------------------------------

! create an ESMF_Array from an AttrVect and a gsMap (for the distgrid)
function mct2esmf_init_Array_AvGSmap(attrVect, gsMap, name, rc)

    implicit none

    ! inout parameters
    type(mct_aVect), intent(in)             :: attrVect
    type(mct_gsMap), intent(in)             :: gsMap
    character(len=*), intent(in), optional  :: name
    integer, intent(out), optional          :: rc

    ! return
    type(ESMF_Array)                        :: mct2esmf_init_Array_AvGSmap

    ! local variables
    type(ESMF_DistGrid)                     :: distgrid
    type(ESMF_ArraySpec)                    :: arrayspec
    integer                                 :: localrc, nfields
    character(len=256)                      :: lname
    character(len=8096)                     :: mct_names
    character(len=*),parameter :: subname = 'mct2esmf_init_Array_AvGSmap'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    lname = 'undefined'
    if (present(name)) then
       lname = trim(name)
    endif

    distgrid = mct2esmf_init_Distgrid_GSmap(gsMap, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    mct_names = mct_aVect_exportRList2c(attrVect)

    mct2esmf_init_Array_AvGSmap =  mct2esmf_init_Array_DistgridList(distgrid,mct_names,name,rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

end function mct2esmf_init_Array_AvGSmap

!--------------------------------------------------------------------------
! create an ESMF_Array from an Array and a gsMap (for the distgrid)
function mct2esmf_init_Array_ArrayGSmap(array, gsMap, name, rc)

    implicit none

    ! inout parameters
    type(ESMF_array), intent(inout)         :: array
    type(mct_gsMap), intent(in)             :: gsMap
    character(len=*), intent(in), optional  :: name
    integer, intent(out), optional          :: rc

    ! return
    type(ESMF_Array)                        :: mct2esmf_init_Array_ArrayGSmap

    ! local variables
    type(ESMF_DistGrid)                     :: distgrid
    type(ESMF_ArraySpec)                    :: arrayspec
    character(len=8096)                     :: mct_names
    integer                                 :: localrc, nfields
    character(len=256)                      :: lname
    character(len=*),parameter :: subname = 'mct2esmf_init_Arrays_ArrayGSmap'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    lname = 'undefined'
    if (present(name)) then
       lname = trim(name)
    endif

    distgrid = mct2esmf_init_Distgrid_GSmap(gsMap, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)
    
    call ESMF_AttributeGet(array, name="mct_names", value=mct_names, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    mct2esmf_init_Array_ArrayGSmap =  mct2esmf_init_Array_DistgridList(distgrid,mct_names,name,rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

end function mct2esmf_init_Array_ArrayGSmap

!--------------------------------------------------------------------------

! copy data from attrVect to ESMF_Array
subroutine mct2esmf_copy_Array_Av(attrVect, array, rc)

    implicit none

    ! inout parameters
    type(mct_aVect), intent(in)             :: attrVect
    type(ESMF_Array), intent(inout)         :: array
    integer, intent(out), optional          :: rc

    ! local variables
    integer                                 :: a_dsize, a_nfields, localrc
    real(ESMF_KIND_R8), pointer             :: a_fptr(:,:)
    integer                                 :: av_dsize, av_nfields
    integer                                 :: i, j, a_off1, a_off2
    character(len=*),parameter :: subname = 'mct2esmf_copy_Array_Av'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array, localDe=0, farrayPtr=a_fptr, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

! tcraig - this call seems to change answers, why!!
!    call esmfshr_util_ArrayGetSize(array, lsize1=a_nfields, lsize2=a_dsize)
    a_dsize = ubound(a_fptr, 2)-lbound(a_fptr, 2)+1
    a_nfields = ubound(a_fptr, 1)-lbound(a_fptr, 1)+1

    av_dsize = mct_aVect_lsize(aV=attrVect) 
    av_nfields = mct_aVect_nRAttr(aV=attrVect)

    if (av_dsize == 0 .or. av_nfields == 0) then
       ! nothing to copy, skip size check too
       return
    endif

    if(a_dsize /= av_dsize .or. a_nfields /= av_nfields) then
        write(shr_log_unit,*) subname,' ERROR: size ',a_dsize,av_dsize,a_nfields,av_nfields
        call shr_sys_flush(shr_log_unit)
        call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)
    endif
        
    ! a_fptr = attrVect%rattr
    a_off2 = lbound(a_fptr, 2)-1
    a_off1 = lbound(a_fptr, 1)-1
    do j = 1, a_dsize
        do i = 1, a_nfields
            a_fptr(a_off1+i, a_off2+j) = attrVect%rattr(i, j)
        enddo
    enddo

end subroutine mct2esmf_copy_Array_Av

!--------------------------------------------------------------------------

! convert Attribute Vectors and Global Segmentation Map to Arrays
function mct2esmf_init_Arrays_AvsGSmap(importAttrVect, exportAttrVect, gsmap, importArray, exportArray, rc)

    implicit none

    ! return
    integer                                 :: mct2esmf_init_Arrays_AvsGSmap

    ! inout parameters
    type(mct_aVect), intent(in)             :: importAttrVect
    type(mct_aVect), intent(in)             :: exportAttrVect
    type(mct_gsMap), intent(in)             :: gsmap
    type(ESMF_Array), intent(inout)         :: importArray
    type(ESMF_Array), intent(inout)         :: exportArray
    integer, intent(out), optional          :: rc

    ! internal variables
    integer                                 :: localrc
    character(len=*),parameter :: subname = 'mct2esmf_init_Arrays_AvsGSmap'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    ! initializations
    importArray = mct2esmf_init_Array_AvGSmap(importAttrVect, gsmap, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    exportArray = mct2esmf_init_Array_AvGSmap(exportAttrVect, gsmap, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    mct2esmf_init_Arrays_AvsGSmap = ESMF_SUCCESS

end function mct2esmf_init_Arrays_AvsGSmap

!--------------------------------------------------------------------------
function mct2esmf_init_Distgrid_Gindex(gindex, rc)

    implicit none

    ! inout parameters
    integer, intent(in)                     :: gindex(:)
    integer, intent(out), optional          :: rc

    ! return
    type(ESMF_DistGrid)                     :: mct2esmf_init_Distgrid_Gindex

    ! local variables
    integer                                 :: localrc
    character(len=*),parameter :: subname = 'mct2esmf_init_Distgrid_Gindex'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    mct2esmf_init_Distgrid_Gindex = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

end function mct2esmf_init_Distgrid_Gindex

!--------------------------------------------------------------------------

! create a 2d Array from index list and colon delimited string
function mct2esmf_init_Array_GindexList(gindex, attname, name, value, rc)

    implicit none

    ! inout parameters
    integer, intent(in)                     :: gindex(:)    ! index list
    character(*), intent(in)                :: attname ! encoded name used later to create mct
    character(*), intent(in), optional      :: name ! encoded name used later to create mct
    real*8, intent(in),optional             :: value  ! initial value
    integer, intent(out),optional           :: rc

    type(ESMF_Array)                        :: mct2esmf_init_Array_Gindexlist

    type(ESMF_DistGrid)                     :: distgrid
    integer                                 :: localrc
    character(len=*),parameter :: subname = 'mct2esmf_init_Array_GindexList'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    ! create 1d arbitrarily distributed distgrid
    distgrid = mct2esmf_init_Distgrid_Gindex(gindex, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    if (         present(name) .and.      present(value)) then
       mct2esmf_init_Array_GindexList =  mct2esmf_init_Array_DistgridList(distgrid,attname,name,value,rc=localrc)
    elseif (     present(name) .and. .not.present(value)) then
       mct2esmf_init_Array_GindexList =  mct2esmf_init_Array_DistgridList(distgrid,attname,name,rc=localrc)
    elseif (.not.present(name) .and.      present(value)) then
       mct2esmf_init_Array_GindexList =  mct2esmf_init_Array_DistgridList(distgrid,attname,value=value,rc=localrc)
    elseif (.not.present(name) .and. .not.present(value)) then
       mct2esmf_init_Array_GindexList =  mct2esmf_init_Array_DistgridList(distgrid,attname,rc=localrc)
    endif

    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

end function mct2esmf_init_Array_GindexList

!--------------------------------------------------------------------------

! create a 2d Array from distgrid and colon delimited string
function mct2esmf_init_Array_DistgridList(distgrid, attname, name, value, rc)

    implicit none

    ! inout parameters
    type(ESMF_DistGrid),intent(in)          :: distgrid ! input decomp
    character(*), intent(in)                :: attname  ! fields list
    character(*), intent(in), optional      :: name     ! name of array
    real*8, intent(in),optional             :: value    ! initial value
    integer, intent(out),optional           :: rc

    type(ESMF_Array)                        :: mct2esmf_init_Array_DistgridList

    ! internal variables
    type(ESMF_ArraySpec)                    :: arrayspec
    real(ESMF_KIND_R8), pointer             :: fptr(:,:)
    integer                                 :: localrc, nfields,lsize1,lsize2
    character(len=256)                      :: lname
    character(len=*),parameter :: subname = 'mct2esmf_init_Array_DistgridList'
    !---------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS

    lname = 'undefined'
    if (present(name)) then
       lname = trim(name)
    endif

    call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=ESMF_TYPEKIND_R8, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    ! convert name to attribute and attach to Array
    
    nfields = shr_string_listGetNum(attname)

    ! create a 2D array, 1d undistributed index of fields, 2d is packed data
    mct2esmf_init_Array_DistgridList = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
        undistLBound=(/1/), undistUBound=(/nfields/), &
        indexflag = ESMF_INDEX_DELOCAL, &
        name=trim(lname), &
        rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    call ESMF_AttributeSet(mct2esmf_init_Array_DistgridList, name="mct_names", value=attname, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    call ESMF_ArrayGet(mct2esmf_init_Array_DistgridList, localDe=0, farrayPtr=fptr, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)
    if (present(value)) then
       fptr = value
    else
       fptr = 0.0_ESMF_KIND_R8
    endif

end function mct2esmf_init_Array_DistgridList

!--------------------------------------------------------------------------

! create an ESMF_Array from an nFields and a gsMap (for the distgrid)
function mct2esmf_init_Array_gsMapNF(nFields, gsMap, name, rc)

    ! inout parameters
    integer, intent(in)                     :: nFields
    type(mct_gsMap), intent(in)             :: gsMap
    character(len=*), intent(in), optional  :: name
    integer, intent(out), optional          :: rc

    ! return
    type(ESMF_Array)                        :: mct2esmf_init_Array_gsMapNF

    ! local variables
    type(ESMF_DistGrid)                     :: distgrid
    type(ESMF_ArraySpec)                    :: arrayspec
    integer                                 :: localrc

    if(present(rc)) rc = ESMF_SUCCESS

    distgrid = mct2esmf_init_DistGrid_GSmap(gsMap, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=ESMF_TYPEKIND_R8, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    ! create a 2D array, 1d undistributed index of fields, 2d is packed data
    mct2esmf_init_Array_gsMapNF = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
        undistLBound=(/1/), undistUBound=(/nfields/), &
        name = name, &
        rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

end function mct2esmf_init_Array_gsMapNF

!--------------------------------------------------------------------------

! create an ESMF_Array from an rList and a gsMap (for the distgrid)
function mct2esmf_init_Array_gsMaprList(gsMap, rList, name, rc)

    ! inout parameters
    type(mct_gsMap), intent(in)             :: gsMap
    character(len=*), intent(in)            :: rList
    character(len=*), intent(in), optional  :: name
    integer, intent(out), optional          :: rc

    ! return
    type(ESMF_Array)                        :: mct2esmf_init_Array_gsMaprList

    ! local variables
    type(ESMF_DistGrid)                     :: distgrid
    type(ESMF_ArraySpec)                    :: arrayspec
    integer                                 :: localrc, nfields
    type(mct_List)                          :: aList

    if(present(rc)) rc = ESMF_SUCCESS

    distgrid = mct2esmf_init_DistGrid_GSmap(gsMap, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    mct2esmf_init_Array_gsMaprList = mct2esmf_init_Array_DistgridList(distgrid, rList, name, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

end function mct2esmf_init_Array_gsMaprList

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------

! copy data from attrVect to ESMF_Array
! src and dst object must have identical distribution, the ArrayCopy routine will check it.
subroutine mct2esmf_copy_alist(attrVect, array, list, rc)

    ! inout parameters
    type(mct_aVect), intent(in)             :: attrVect
    type(ESMF_Array), intent(inout)         :: array
    character(len=*), intent(in)            :: list
    integer, intent(out), optional          :: rc

    ! local variables
    type(ESMF_Array)                        :: srcArray
    type(ESMF_DistGrid)                     :: distgrid
    character(len=4096)                     :: flist
    integer                                 :: localrc

    if(present(rc)) rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array, distgrid=distgrid, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    !--------------------------------------------------------------------------
    ! use the attribute vector Fortran array pointer to create a temporary Array;
    ! this avoids unnecessary memory copy and code is more efficient
    !--------------------------------------------------------------------------
    srcArray = ESMF_ArrayCreate(farrayPtr=attrVect%rattr, distgrid=distgrid, distgridToArrayMap=(/2/), &
        rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    call ESMF_AttributeSet(srcArray, name="mct_names", &
        value=trim(mct_aVect_exportRList2c(attrVect)), rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    ! do a field name based copy between incongruent src and dst Array pairs
    call esmfshr_util_ArrayCopy(srcArray, array, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

    call ESMF_ArrayDestroy(srcArray, rc=localrc)
    if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, terminationflag=ESMF_ABORT)

end subroutine mct2esmf_copy_alist

#endif
end module
