!-------------------------------------------------------------------
! rebins the 4 sea salt bins into 2 bins for the radiation
!
!  N.B. This code looks for the constituents of SSLTA and SSLTC
!       in the physics buffer first, and uses those if found.
!       Consequently, it is not possible to have prognostic sea
!       salt be radiatively active if the prescribed sea salt is
!       also present.  The current (cam3_5_52) chemistry configurations
!       don't allow both prescribed and prognostic to be present
!       simultaneously, but a more flexible chemistry package that
!       allows this would break this code.
!
! Created by: Francis Vitt
! Date: 9 May 2008
!-------------------------------------------------------------------
module sslt_rebin

  use shr_kind_mod,   only: r8 => shr_kind_r8

  implicit none

  integer :: indices(4)
  integer :: sslta_idx, ssltc_idx

  logical :: has_sslt = .false.
  character(len=1) :: source
  character(len=1), parameter :: DATA = 'D'
  character(len=1), parameter :: PROG = 'P'

  private
  public :: sslt_rebin_init, sslt_rebin_adv, sslt_rebin_register
contains


!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sslt_rebin_register
    use ppgrid,       only : pver
    use phys_buffer,  only : pbuf_add

    ! add SSLTA and SSLTC to physics buffer
    call pbuf_add('SSLTA','physpkg',1,pver,1,sslta_idx)
    call pbuf_add('SSLTC','physpkg',1,pver,1,ssltc_idx)

  endsubroutine sslt_rebin_register

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sslt_rebin_init

    use constituents, only : cnst_get_ind
    use phys_buffer,  only : pbuf, pbuf_get_fld_idx
    use ppgrid,       only : pver
    use cam_history,  only : addfld, phys_decomp

    implicit none

    indices(1) = pbuf_get_fld_idx('sslt1', failcode=-1 )
    indices(2) = pbuf_get_fld_idx('sslt2', failcode=-1 )
    indices(3) = pbuf_get_fld_idx('sslt3', failcode=-1 )
    indices(4) = pbuf_get_fld_idx('sslt4', failcode=-1 )

    has_sslt = all( indices(:) > 0 )
    if ( has_sslt ) source = DATA

    if ( .not. has_sslt ) then
       call cnst_get_ind ('SSLT01', indices(1), abort=.false.)
       call cnst_get_ind ('SSLT02', indices(2), abort=.false.)
       call cnst_get_ind ('SSLT03', indices(3), abort=.false.)
       call cnst_get_ind ('SSLT04', indices(4), abort=.false.)
       has_sslt = all( indices(:) > 0 )
       if ( has_sslt ) source = PROG
    endif

    if ( has_sslt ) then
       call addfld('SSLTA','kg/kg', pver, 'A', 'sea salt', phys_decomp )
       call addfld('SSLTC','kg/kg', pver, 'A', 'sea salt', phys_decomp )
    endif

    ! initialize the pbuf values
    pbuf(sslta_idx)%fld_ptr = 0._r8
    pbuf(ssltc_idx)%fld_ptr = 0._r8

  end subroutine sslt_rebin_init
  
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sslt_rebin_adv(pbuf, phys_state)

    use physics_types,only : physics_state
    use phys_buffer,  only : pbuf_size_max, pbuf_fld
    use ppgrid,       only : pver, pcols
    use cam_history,  only : outfld

    implicit none

    type(pbuf_fld),      intent(in) :: pbuf(pbuf_size_max)
    type(physics_state), target, intent(in) :: phys_state

!++ changed wgt_sscm declaration for roundoff validation with earlier code
!    real(r8), parameter :: wgt_sscm = 6.0_r8 / 7.0_r8 ! Fraction of total seasalt mass in coarse mode 
    real, parameter :: wgt_sscm = 6.0 / 7.0 ! Fraction of total seasalt mass in coarse mode 

    real(r8), dimension(:,:), pointer :: sslt1, sslt2, sslt3, sslt4
    real(r8), dimension(:,:), pointer :: sslta, ssltc
    integer :: lchnk, ncol
    real(r8) :: sslt_sum(pcols,pver)

    lchnk = phys_state%lchnk
    ncol = phys_state%ncol

    if (.not. has_sslt) return

    select case( source )
    case (PROG)
       sslt1 => phys_state%q(:,:,indices(1))
       sslt2 => phys_state%q(:,:,indices(2))
       sslt3 => phys_state%q(:,:,indices(3))
       sslt4 => phys_state%q(:,:,indices(4))
    case (DATA)
       sslt1 => pbuf(indices(1))%fld_ptr(1,:,:,lchnk,1)
       sslt2 => pbuf(indices(2))%fld_ptr(1,:,:,lchnk,1)
       sslt3 => pbuf(indices(3))%fld_ptr(1,:,:,lchnk,1)
       sslt4 => pbuf(indices(4))%fld_ptr(1,:,:,lchnk,1)
    end select

    sslta => pbuf(sslta_idx)%fld_ptr(1,:,:,lchnk,1)
    ssltc => pbuf(ssltc_idx)%fld_ptr(1,:,:,lchnk,1)

    sslt_sum(:ncol,:) = sslt1(:ncol,:) + sslt2(:ncol,:) + sslt3(:ncol,:) + sslt4(:ncol,:)
    sslta(:ncol,:) = (1._r8-wgt_sscm)*sslt_sum(:ncol,:) ! fraction of seasalt mass in accumulation mode
    ssltc(:ncol,:) = wgt_sscm*sslt_sum(:ncol,:) ! fraction of seasalt mass in coagulation mode

    call outfld( 'SSLTA', sslta(:ncol,:), ncol, lchnk )
    call outfld( 'SSLTC', ssltc(:ncol,:), ncol, lchnk )

  end subroutine sslt_rebin_adv

end module sslt_rebin
