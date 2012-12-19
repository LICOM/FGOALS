module map_atmice_mct

!---------------------------------------------------------------------
!
! Purpose:
!
! Collect coupling routines for sequential coupling of ICE-ATM.
!
! Author: R. Jacob, M. Vertenstein
!
!---------------------------------------------------------------------

  use shr_kind_mod     , only: R8 => SHR_KIND_R8
  use shr_sys_mod

  use shr_mct_mod, only: shr_mct_sMatPInitnc
  use mct_mod
  use seq_flds_mod
  use seq_flds_indices
  use seq_cdata_mod
  use seq_comm_mct, only: logunit, loglevel
  use seq_infodata_mod
  use seq_map_mod
  use m_die

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: map_ice2atm_init_mct
  public :: map_ice2atm_mct
! atm2ice is not used or validated yet
  private :: map_atm2ice_init_mct   
  private :: map_atm2ice_mct

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(mct_rearr), private :: Re_ice2atm
  type(mct_sMatp), private :: sMatp_Fi2a
  type(mct_sMatp), private :: sMatp_Si2a
! atm2ice is not used or validated yet
  type(mct_rearr), private :: Re_atm2ice
  type(mct_sMatp), private :: sMatp_Fa2i
  type(mct_sMatp), private :: sMatp_Sa2i

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
  logical :: samegrid_mapa2i

!=======================================================================
contains
!=======================================================================

  subroutine map_atm2ice_init_mct( cdata_a, cdata_i)

    !-----------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_a
    type(seq_cdata),intent(in) :: cdata_i
    !
    ! Local variables
    !
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap), pointer :: gsMap_a
    type(mct_gsMap), pointer :: gsMap_i
    integer                  :: mpicom
    character(*),parameter :: subName = '(map_atm2ice_init_mct) '
    !-----------------------------------------------

    call seq_cdata_setptrs(cdata_a, gsMap=gsMap_a)
    call seq_cdata_setptrs(cdata_i, gsMap=gsMap_i)
    call seq_cdata_setptrs(cdata_a, mpicom=mpicom, infodata=infodata)

    call seq_infodata_GetData( infodata, samegrid_ao=samegrid_mapa2i)

    if (samegrid_mapa2i) then

       call mct_rearr_init(gsMap_a, gsMap_i, mpicom, Re_atm2ice)

    else

       call shr_mct_sMatPInitnc(sMatp_Fa2i,gsMap_a,gsMap_i,"seq_maps.rc", &
            "atm2iceFmapname:","atm2iceFmaptype:",mpicom)
       call shr_mct_sMatPInitnc(sMatp_Sa2i,gsMap_a,gsMap_i,"seq_maps.rc", &
            "atm2iceSmapname:","atm2iceSmaptype:",mpicom)

    endif
    
  end subroutine map_atm2ice_init_mct

!=======================================================================

  subroutine map_ice2atm_init_mct( cdata_i, cdata_a)

    !-----------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in) :: cdata_i
    type(seq_cdata), intent(in) :: cdata_a
    !
    ! Local variables
    !
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap), pointer :: gsMap_i           ! ice gsMap
    type(mct_gsMap), pointer :: gsMap_a           ! atm gsMap
    type(mct_gGrid), pointer :: dom_i             ! ice domain
    type(mct_gGrid), pointer :: dom_a             ! atm domain
    integer                  :: mpicom            ! communicator spanning atm and ice
    character(*),parameter :: subName = '(map_ice2atm_init_mct) '
    !-----------------------------------------------

    call seq_cdata_setptrs(cdata_i, gsMap=gsMap_i, dom=dom_i)
    call seq_cdata_setptrs(cdata_a, gsMap=gsMap_a, dom=dom_a)
    call seq_cdata_setptrs(cdata_a, mpicom=mpicom, infodata=infodata)

    call seq_infodata_GetData( infodata, samegrid_ao=samegrid_mapa2i)

    ! Initialize ice->atm mapping or rearranging

    if(samegrid_mapa2i) then

      call mct_rearr_init(gsMap_i, gsMap_a, mpicom, Re_ice2atm)

    else

      call shr_mct_sMatPInitnc(sMatp_Fi2a, gsMap_i, gsMap_a, "seq_maps.rc", &
        "ice2atmFmapname:", "ice2atmFmaptype:", mpicom)

      call shr_mct_sMatPInitnc(sMatp_Si2a, gsMap_i, gsMap_a, "seq_maps.rc", &
        "ice2atmSmapname:", "ice2atmSmaptype:", mpicom)

    endif

  end subroutine map_ice2atm_init_mct

!=======================================================================

  subroutine map_atm2ice_mct( cdata_a, av_a, cdata_i, av_i, &
                              fluxlist, statelist)

    !-----------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata) ,intent(in)           :: cdata_a
    type(mct_aVect) ,intent(in)           :: av_a
    type(seq_cdata) ,intent(in)           :: cdata_i
    type(mct_aVect) ,intent(out)          :: av_i
    character(len=*),intent(in), optional :: fluxlist
    character(len=*),intent(in), optional :: statelist
    !
    ! Local variables
    !
    integer         :: lsize      ! temporary
    type(mct_aVect) :: av_i_f     ! temporary flux attribute vector
    type(mct_aVect) :: av_i_s     ! temporary state attribute vector
    type(mct_aVect) :: av_a_fl    ! temporary av for rearranging
    type(mct_aVect) :: av_i_fl    ! temporary av for rearranging
    character(*),parameter :: subName = '(map_atm2ice_mct) '
    !-----------------------------------------------

    if(samegrid_mapa2i) then

      if (present(fluxlist) .or. present(statelist)) then
         if (present(fluxlist)) then
            call mct_rearr_rearrange_fldlist(av_a, av_i, Re_atm2ice, VECTOR=usevector, &
                ALLTOALL=usealltoall, fldlist=fluxlist)
         end if
         if (present(statelist)) then
            call mct_rearr_rearrange_fldlist(av_a, av_i, Re_atm2ice, VECTOR=usevector, &
                 ALLTOALL=usealltoall, fldlist=statelist)
         end if
      else
         call mct_rearr_rearrange(av_a, av_i, Re_atm2ice, VECTOR=usevector, ALLTOALL=usealltoall)
      endif

    else

      if (present(fluxlist) .or. present(statelist)) then
         if (present(fluxlist)) then
!            lsize = mct_aVect_lsize(av_i)
!            call mct_aVect_init (av_i_f, rlist=fluxlist , lsize=lsize)
!            call mct_sMat_avMult(av_a, sMatp_Fa2i, av_i_f, VECTOR=usevector, rList=fluxlist)
!            call mct_aVect_copy (aVin=av_i_f, aVout=av_i, vector=usevector)
!            call mct_aVect_clean(av_i_f)
            call seq_map_avNorm(av_a, av_i, sMatp_Fa2i, rList=fluxlist, donorm=.false.)
         end if
         if (present(statelist)) then
!            lsize = mct_aVect_lsize(av_i)
!            call mct_aVect_init (av_i_s, rlist=statelist , lsize=lsize)
!            call mct_sMat_avMult(av_a, sMatp_Sa2i, av_i_s, VECTOR=usevector, rList=statelist)
!            call mct_aVect_copy (aVin=av_i_s, aVout=av_i, vector=usevector)
!            call mct_aVect_clean(av_i_s)
            call seq_map_avNorm(av_a, av_i, sMatp_Sa2i, rList=statelist, donorm=.false.)
         end if
      else
         !--- default is flux mapping
!         call mct_sMat_avMult(av_a, sMatp_Fa2i, av_i, VECTOR=usevector)
         call seq_map_avNorm(av_a, av_i, sMatp_Fa2i, donorm=.false.)
      endif

    endif

  end subroutine map_atm2ice_mct

!=======================================================================

  subroutine map_ice2atm_mct( cdata_i, av_i, cdata_a, av_a, &
	                      fractions_i, fractions_a, &
                              fluxlist, statelist)

    !-----------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata) ,intent(in)           :: cdata_i
    type(mct_AVect) ,intent(in)           :: av_i
    type(seq_cdata) ,intent(in)           :: cdata_a
    type(mct_AVect) ,intent(out)          :: av_a
    type(mct_AVect) ,intent(in), optional :: fractions_i
    type(mct_AVect) ,intent(in), optional :: fractions_a
    character(len=*),intent(in), optional :: fluxlist
    character(len=*),intent(in), optional :: statelist
    !
    ! Local Variables
    !
    integer              :: n,ki,i,j       ! indices
    integer              :: numats,ier     ! number of attributes
    integer              :: lsize          ! size of attribute vector
    type(mct_aVect)      :: av_a_f, av_a_s ! temporary
    character(*),parameter :: subName = '(map_ice2atm_mct) '
    !-----------------------------------------------

    if (samegrid_mapa2i) then

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             call mct_rearr_rearrange_fldlist(av_i, av_a, Re_ice2atm, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          end if
          if (present(statelist)) then
             call mct_rearr_rearrange_fldlist(av_i, av_a, Re_ice2atm, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          end if
       else
          call mct_rearr_rearrange(av_i, av_a, Re_ice2atm, VECTOR=usevector, ALLTOALL=usealltoall)
       endif

    else

       if (present(fractions_i)) then

          ! Perform mapping

          if (present(fluxlist) .or. present(statelist)) then
             if (present(fluxlist)) then 
                call seq_map_avNorm(av_i, av_a, sMatp_Fi2a, fractions_i, 'ifrac', fractions_a, 'ifrac', rList=fluxlist)
             end if
             if (present(statelist)) then
                call seq_map_avNorm(av_i, av_a, sMatp_Si2a, fractions_i, 'ifrac', fractions_a, 'ifrac', rList=statelist)
             end if
          else
             call seq_map_avNorm(av_i, av_a, sMatp_Fi2a, fractions_i, 'ifrac', fractions_a, 'ifrac')
          endif
          
       else

          if (present(fluxlist) .or. present(statelist)) then
             if (present(fluxlist)) then 
!                lsize = mct_aVect_lsize(av_a)
!                call mct_aVect_init (av_a_f, rlist=fluxlist , lsize=lsize)
!                call mct_sMat_avMult(av_i, sMatp_Fi2a, av_a_f, VECTOR=usevector, rList=fluxlist)
!                call mct_aVect_copy (aVin=av_a_f, aVout=av_a, vector=usevector)
!                call mct_aVect_clean(av_a_f)
                call seq_map_avNorm(av_i, av_a, sMatp_Fi2a, rList=fluxlist, donorm=.false.)
             end if
             if (present(statelist)) then
!                lsize = mct_aVect_lsize(av_a)
!                call mct_aVect_init (av_a_s, rlist=statelist, lsize=lsize)
!                call mct_sMat_avMult(av_i, sMatp_Si2a, av_a_s, VECTOR=usevector, rList=statelist)
!                call mct_aVect_copy (aVin=av_a_s, aVout=av_a, vector=usevector)
!                call mct_aVect_clean(av_a_s)
                call seq_map_avNorm(av_i, av_a, sMatp_Si2a, rList=statelist, donorm=.false.)
             end if
          else
             !--- default is flux mapping
!             call mct_sMat_avMult(av_i, sMatp_Fi2a, av_a, VECTOR=usevector)
             call seq_map_avNorm(av_i, av_a, sMatp_Fi2a, donorm=.false.)
          endif
         
      end if

   endif

 end subroutine map_ice2atm_mct

!=======================================================================

end module map_atmice_mct
