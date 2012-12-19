module cloud_diagnostics

!---------------------------------------------------------------------------------
! Purpose:
!
! Put cloud physical specifications on the history tape
!  Modified from code that computed cloud optics
!
! Author: Byron Boville  Sept 06, 2002
!  Modified Oct 15, 2008
!    
!
!---------------------------------------------------------------------------------

   use shr_kind_mod,  only: r8=>shr_kind_r8
   use ppgrid,        only: pcols, pver
   use physconst,     only: gravit

   implicit none
   private
   save

   public :: cloud_diagnostics_init, put_cloud_diagnostics

! Local variables
   integer :: &
     i_dei, i_mu, i_lambda, i_iciwp, i_iclwp, i_cld  ! index into pbuf for cloud fields

contains

!===============================================================================
  subroutine cloud_diagnostics_init()
!-----------------------------------------------------------------------
    use cam_history,   only: addfld, add_default, phys_decomp
    use phys_buffer,   only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx, pbuf_times

    implicit none
!-----------------------------------------------------------------------

! register history variables

    call addfld ('ICLDIWP', 'kg/m2', pver, 'A','In-cloud ice water path'               ,phys_decomp)
    call addfld ('ICLDTWP ','kg/m2 ',pver, 'A','In-cloud cloud total water path (liquid and ice)',phys_decomp)

    call addfld ('GCLDLWP ','kg/m2 ',pver, 'A','Grid-box cloud water path'             ,phys_decomp)
    call addfld ('TGCLDCWP','kg/m2 ',1,    'A','Total grid-box cloud water path (liquid and ice)',phys_decomp)
    call addfld ('TGCLDLWP','kg/m2 ',1,    'A','Total grid-box cloud liquid water path',phys_decomp)
    call addfld ('TGCLDIWP','kg/m2 ',1,    'A','Total grid-box cloud ice water path'   ,phys_decomp)
    call add_default ('TGCLDLWP', 1, ' ')
    call add_default ('TGCLDIWP', 1, ' ')
    call add_default ('TGCLDCWP', 1, ' ')

    call addfld ('CLOUD   ','fraction',pver, 'A','Cloud fraction'                        ,phys_decomp)
    call addfld ('CLDTOT  ','fraction',1,    'A','Vertically-integrated total cloud'     ,phys_decomp)
    call addfld ('CLDLOW  ','fraction',1,    'A','Vertically-integrated low cloud'       ,phys_decomp)
    call addfld ('CLDMED  ','fraction',1,    'A','Vertically-integrated mid-level cloud' ,phys_decomp)
    call addfld ('CLDHGH  ','fraction',1,    'A','Vertically-integrated high cloud'      ,phys_decomp)
    call add_default ('CLOUD   ', 1, ' ')
    call add_default ('CLDTOT  ', 1, ' ')
    call add_default ('CLDLOW  ', 1, ' ')
    call add_default ('CLDMED  ', 1, ' ')
    call add_default ('CLDHGH  ', 1, ' ')

    i_dei    = pbuf_get_fld_idx('DEI')
    i_mu     = pbuf_get_fld_idx('MU')
    i_lambda = pbuf_get_fld_idx('LAMBDAC')

    i_iciwp  = pbuf_get_fld_idx('ICIWP')
    i_iclwp  = pbuf_get_fld_idx('ICLWP')
    i_cld    = pbuf_get_fld_idx('CLD')

    call addfld ('lambda_cloud','1/meter',pver,'I','lambda in cloud', phys_decomp)
    call addfld ('mu_cloud','1',pver,'I','mu in cloud', phys_decomp)
    call addfld ('dei_cloud','micrometers',pver,'I','ice radiative effective diameter in cloud', phys_decomp)

    call addfld ('SETLWP  ','gram/m2 ',pver, 'A','Prescribed liquid water path'          ,phys_decomp)
    call addfld ('LWSH    ','m       ',1,    'A','Liquid water scale height'             ,phys_decomp)

    return
  end subroutine cloud_diagnostics_init

!===============================================================================
  subroutine put_cloud_diagnostics(state, pbuf)
!
! Compute (liquid+ice) water path and cloud water/ice diagnostics
! *** soon this code will compute liquid and ice paths from input liquid and ice mixing ratios
! 
! **** mixes interface and physics code temporarily
!-----------------------------------------------------------------------
    use physics_types, only: physics_state
    use phys_buffer,   only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx, pbuf_times
    use cam_history,   only: outfld
    use pkg_cldoptics, only: cldovrlap, cldclw

    implicit none

! Arguments
    type(physics_state), intent(in)  :: state        ! state variables
    type(pbuf_fld),      intent(in), dimension(pbuf_size_max) :: pbuf

    real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction
    real(r8), pointer, dimension(:,:) :: iciwpth    ! in-cloud cloud ice water path
    real(r8), pointer, dimension(:,:) :: iclwpth    ! in-cloud cloud liquid water path
    real(r8), pointer, dimension(:,:) :: dei        ! effective radiative diameter of ice
    real(r8), pointer, dimension(:,:) :: mu         ! gamma distribution for liq clouds
    real(r8), pointer, dimension(:,:) :: lambda     ! gamma distribution for liq clouds
    integer :: itim
! Local variables
    real(r8) :: pmxrgn(pcols,pver+1)    ! Maximum values of pressure for each
    integer  :: nmxrgn(pcols)           ! Number of maximally overlapped regions

    real(r8) :: cwp   (pcols,pver)                   ! in-cloud cloud (total) water path
    real(r8) :: cicewp(pcols,pver)                   ! in-cloud cloud ice water path
    real(r8) :: cliqwp(pcols,pver)                   ! in-cloud cloud liquid water path
    real(r8) :: gicewp(pcols,pver)                   ! grid-box cloud ice water path
    real(r8) :: gliqwp(pcols,pver)                   ! grid-box cloud liquid water path
    real(r8) :: gwp   (pcols,pver)                   ! grid-box cloud (total) water path
    real(r8) :: tgicewp(pcols)                       ! Vertically integrated ice water path
    real(r8) :: tgliqwp(pcols)                       ! Vertically integrated liquid water path
    real(r8) :: tgwp   (pcols)                       ! Vertically integrated (total) cloud water path

! old data
    real(r8) :: tpw    (pcols)                       ! total precipitable water
    real(r8) :: clwpold(pcols,pver)                  ! Presribed cloud liq. h2o path
    real(r8) :: hl     (pcols)                       ! Liquid water scale height

    real(r8) cltot(pcols)                      ! Diagnostic total cloud cover
    real(r8) cllow(pcols)                      !       "     low  cloud cover
    real(r8) clmed(pcols)                      !       "     mid  cloud cover
    real(r8) clhgh(pcols)                      !       "     hgh  cloud cover

    integer :: i,k                                   ! loop indexes
    integer :: ncol, lchnk
    real(r8) :: rgrav
!-----------------------------------------------------------------------
    ncol  = state%ncol
    lchnk = state%lchnk

    itim = pbuf_old_tim_idx()
    cld     => pbuf(i_cld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
    iclwpth => pbuf(i_iclwp)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
    iciwpth => pbuf(i_iciwp)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
    dei     => pbuf(i_dei)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
    mu      => pbuf(i_mu)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
    lambda  => pbuf(i_lambda)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

! Compute liquid and ice water paths
    tgicewp(:ncol) = 0._r8
    tgliqwp(:ncol) = 0._r8
    do k=1,pver
       do i = 1,ncol
          gicewp(i,k) = iciwpth(i,k)*cld(i,k)
          gliqwp(i,k) = iclwpth(i,k)*cld(i,k)
          cicewp(i,k) = iciwpth(i,k)
          cliqwp(i,k) = iclwpth(i,k)
          tgicewp(i)  = tgicewp(i) + gicewp(i,k)
          tgliqwp(i)  = tgliqwp(i) + gliqwp(i,k)
       end do
    end do
    tgwp(:ncol) = tgicewp(:ncol) + tgliqwp(:ncol)
    gwp(:ncol,:pver) = gicewp(:ncol,:pver) + gliqwp(:ncol,:pver) 
    cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver) 

! Compute total preciptable water in column (in mm)
    tpw(:ncol) = 0.0_r8
    rgrav = 1.0_r8/gravit
    do k=1,pver
       do i=1,ncol
          tpw(i) = tpw(i) + state%pdel(i,k)*state%q(i,k,1)*rgrav
       end do
    end do

! Determine parameters for maximum/random overlap
    call cldovrlap(lchnk, ncol, state%pint, cld, nmxrgn, pmxrgn)

    call cldsav (lchnk, ncol, cld, state%pmid, cltot, &
            cllow, clmed, clhgh, nmxrgn, pmxrgn)
    !
    ! Dump cloud field information to history tape buffer (diagnostics)
    !
    call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
    call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
    call outfld('CLDMED  ',clmed  ,pcols,lchnk)
    call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)

    call outfld('CLOUD   ',cld    ,pcols,lchnk) 

    call outfld('GCLDLWP' ,gwp    , pcols,lchnk)
    call outfld('TGCLDCWP',tgwp   , pcols,lchnk)
    call outfld('TGCLDLWP',tgliqwp, pcols,lchnk)
    call outfld('TGCLDIWP',tgicewp, pcols,lchnk)
    call outfld('ICLDTWP' ,cwp    , pcols,lchnk)
    call outfld('ICLDIWP' ,cicewp , pcols,lchnk)

    call outfld('dei_cloud',dei(:,:),pcols,lchnk)
    call outfld('mu_cloud',mu(:,:),pcols,lchnk)
    call outfld('lambda_cloud',lambda(:,:),pcols,lchnk)

! Diagnostic liquid water path (old specified form)

    call cldclw(lchnk, ncol, state%zi, clwpold, tpw, hl)
    call outfld('SETLWP'  ,clwpold, pcols,lchnk)
    call outfld('LWSH'    ,hl     , pcols,lchnk)

    return
end subroutine put_cloud_diagnostics


end module cloud_diagnostics
