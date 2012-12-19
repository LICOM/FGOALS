!================================================================================================
! output data necessary to drive radiation offline
! Francis Vitt -- Created 15 Dec 2009
!================================================================================================
module radiation_data

  use shr_kind_mod,     only: r8=>shr_kind_r8
  use ppgrid,           only: pcols, pver, pverp
  use cam_history,      only: addfld, add_default, phys_decomp, outfld
  use phys_buffer,      only: pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx, pbuf_size_max
  use cam_history_support, only: fieldname_len, fillvalue
  use spmd_utils,       only: masterproc
  use abortutils,       only: endrun

  implicit none
  private

  public :: output_rad_data
  public :: init_rad_data
  public :: rad_data_readnl

  integer, public :: cld_ifld,concld_ifld,rel_ifld,rei_ifld
  integer, public :: dei_ifld,mu_ifld,lambdac_ifld,iciwp_ifld,iclwp_ifld,rel_fn_ifld
  integer, public :: des_ifld,icswp_ifld,cldfsnow_ifld

  character(len=fieldname_len), public, parameter :: &
       lndfrc_fldn    = 'rad_lndfrc      ' , &
       icefrc_fldn    = 'rad_icefrc      ' , &
       snowh_fldn     = 'rad_snowh       ' , &
       landm_fldn     = 'rad_landm       ' , &
       asdir_fldn     = 'rad_asdir       ' , &
       asdif_fldn     = 'rad_asdif       ' , &
       aldir_fldn     = 'rad_aldir       ' , &
       aldif_fldn     = 'rad_aldif       ' , &
       coszen_fldn    = 'rad_coszen      ' , &
       asdir_pos_fldn = 'rad_asdir_pos   ' , &
       asdif_pos_fldn = 'rad_asdif_pos   ' , &
       aldir_pos_fldn = 'rad_aldir_pos   ' , &
       aldif_pos_fldn = 'rad_aldif_pos   ' , &
       lwup_fldn      = 'rad_lwup        ' , &
       ts_fldn        = 'rad_ts          ' , &
       temp_fldn      = 'rad_temp        ' , &
       pdel_fldn      = 'rad_pdel        ' , &
       pdeldry_fldn   = 'rad_pdeldry     ' , &
       pmid_fldn      = 'rad_pmid        ' , &
       watice_fldn    = 'rad_watice      ' , &
       watliq_fldn    = 'rad_watliq      ' , &
       watvap_fldn    = 'rad_watvap      ' , &
       zint_fldn      = 'rad_zint        ' , &
       pint_fldn      = 'rad_pint        ' , &
       cld_fldn       = 'rad_cld         ' , &
       cldfsnow_fldn  = 'rad_cldfsnow    ' , &
       concld_fldn    = 'rad_concld      ' , &
       rel_fldn       = 'rad_rel         ' , &
       rei_fldn       = 'rad_rei         ' , &
       dei_fldn       = 'rad_dei         ' , &
       des_fldn       = 'rad_des         ' , &
       mu_fldn        = 'rad_mu          ' , &
       lambdac_fldn   = 'rad_lambdac     ' , &
       iciwp_fldn     = 'rad_iciwp       ' , &
       iclwp_fldn     = 'rad_iclwp       ' , &
       icswp_fldn     = 'rad_icswp       ' , &
       rel_fn_fldn    = 'rad_rel_fn      ' 

  ! rad constituents mixing ratios
  integer,           public :: nrad_cnsts
  character(len=64), public, allocatable :: rad_cnstnames(:)
  character(len=1 ), public, allocatable :: rad_cnstsources(:)
  integer,           public, allocatable :: rad_cnstindices(:)
 
  ! control options  
  logical          :: rad_data_output = .false.
  integer          :: rad_data_histfile_num = 2
  character(len=1) :: rad_data_avgflag = 'A'

  ! MG microphys check
  logical, public :: mg_microphys

contains

!================================================================================================
!================================================================================================
  subroutine rad_data_readnl(nlfile)

    ! Read rad_data_nl namelist group.  Parse input.

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr, i
    character(len=*), parameter :: subname = 'rad_data_readnl'

    namelist /rad_data_nl/ rad_data_output, rad_data_histfile_num, rad_data_avgflag

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'rad_data_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, rad_data_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast (rad_data_output,       1,   mpilog ,  0, mpicom)
    call mpibcast (rad_data_histfile_num, 1,   mpiint ,  0, mpicom)
    call mpibcast (rad_data_avgflag,      1,   mpichar , 0, mpicom)
#endif
    
  end subroutine rad_data_readnl

  !================================================================================================
  !================================================================================================
  subroutine init_rad_data
    use rad_constituents, only: rad_cnst_get_info
    use phys_control,     only: phys_getopts

    implicit none
    
    integer :: i, naer, ngas
    character(len=64) :: name
    character(len=128):: long_name
    character(len=64) :: long_name_description
    character(len=16)  :: microp_scheme  ! microphysics scheme

    if (.not.rad_data_output) return
   
    call phys_getopts(microp_scheme_out=microp_scheme)
    mg_microphys =  trim(microp_scheme) .eq. 'MG' 

    cld_ifld     = pbuf_get_fld_idx('CLD')
    concld_ifld  = pbuf_get_fld_idx('CONCLD')
    rel_ifld     = pbuf_get_fld_idx('REL')
    rei_ifld     = pbuf_get_fld_idx('REI')
    if (mg_microphys) then
       dei_ifld     = pbuf_get_fld_idx('DEI')
       des_ifld     = pbuf_get_fld_idx('DES')
       mu_ifld      = pbuf_get_fld_idx('MU')
       lambdac_ifld = pbuf_get_fld_idx('LAMBDAC')
       iciwp_ifld   = pbuf_get_fld_idx('ICIWP')
       iclwp_ifld   = pbuf_get_fld_idx('ICLWP')
       icswp_ifld   = pbuf_get_fld_idx('ICSWP')
       rel_fn_ifld  = pbuf_get_fld_idx('REL_FN')
       cldfsnow_ifld= pbuf_get_fld_idx('CLDFSNOW')
    endif

    call addfld (lndfrc_fldn, 'fraction', 1,    rad_data_avgflag,&
         'radiation input: land fraction',phys_decomp)
    call addfld (icefrc_fldn, 'fraction', 1,    rad_data_avgflag,&
         'radiation input: ice fraction',phys_decomp)
    call addfld (snowh_fldn,  'm',        1,    rad_data_avgflag,&
         'radiation input: water equivalent snow depth',phys_decomp)
    call addfld (landm_fldn,  'none',     1,    rad_data_avgflag,&
         'radiation input: land mask: ocean(0), continent(1), transition(0-1)',phys_decomp)

    call addfld (asdir_fldn,  '1',        1,    rad_data_avgflag,&
         'radiation input: short wave direct albedo',phys_decomp, flag_xyfill=.true.)
    call addfld (asdif_fldn,  '1',        1,    rad_data_avgflag,&
         'radiation input: short wave difuse albedo',phys_decomp, flag_xyfill=.true.)
    call addfld (aldir_fldn,  '1',        1,    rad_data_avgflag,&
         'radiation input: long wave direct albedo', phys_decomp, flag_xyfill=.true.)
    call addfld (aldif_fldn,  '1',        1,    rad_data_avgflag,&
         'radiation input: long wave difuse albedo', phys_decomp, flag_xyfill=.true.)

    call addfld (coszen_fldn,     '1', 1,    rad_data_avgflag,&
         'radiation input: cosine solar zenith when positive', phys_decomp, flag_xyfill=.true.)
    call addfld (asdir_pos_fldn,  '1', 1,    rad_data_avgflag,&
         'radiation input: short wave direct albedo weighted by coszen', phys_decomp, flag_xyfill=.true.)
    call addfld (asdif_pos_fldn,  '1', 1,    rad_data_avgflag,&
         'radiation input: short wave difuse albedo weighted by coszen', phys_decomp, flag_xyfill=.true.)
    call addfld (aldir_pos_fldn,  '1', 1,    rad_data_avgflag,&
         'radiation input: long wave direct albedo weighted by coszen', phys_decomp, flag_xyfill=.true.)
    call addfld (aldif_pos_fldn,  '1', 1,    rad_data_avgflag,&
         'radiation input: long wave difuse albedo weighted by coszen', phys_decomp, flag_xyfill=.true.)
    
    call addfld (lwup_fldn,   'W/m2',     1,    rad_data_avgflag,&
         'radiation input: long wave up radiation flux ',phys_decomp)
    call addfld (ts_fldn,     'K',        1,    rad_data_avgflag,&
         'radiation input: surface temperature',phys_decomp)

    call addfld (temp_fldn,   'K',        pver, rad_data_avgflag,&
         'radiation input: midpoint temperature',phys_decomp)
    call addfld (pdel_fldn,   'Pa',       pver, rad_data_avgflag,&
         'radiation input: pressure layer thickness',phys_decomp)
    call addfld (pdeldry_fldn,'Pa',       pver, rad_data_avgflag,&
         'radiation input: dry pressure layer thickness',phys_decomp)
    call addfld (pmid_fldn,   'Pa',       pver, rad_data_avgflag,&
         'radiation input: midpoint pressure',phys_decomp)
    call addfld (watice_fldn, 'kg/kg',    pver, rad_data_avgflag,&
         'radiation input: cloud ice',phys_decomp)
    call addfld (watliq_fldn, 'kg/kg',    pver, rad_data_avgflag,&
         'radiation input: cloud liquid water',phys_decomp)
    call addfld (watvap_fldn, 'kg/kg',    pver, rad_data_avgflag,&
         'radiation input: water vapor',phys_decomp)

    call addfld (zint_fldn,   'km',       pverp,rad_data_avgflag,&
         'radiation input: interface height',phys_decomp)
    call addfld (pint_fldn,   'Pa',       pverp,rad_data_avgflag,&
         'radiation input: interface pressure',phys_decomp)

    call addfld (cld_fldn,    'fraction', pver, rad_data_avgflag,&
         'radiation input: cloud fraction',phys_decomp)
    call addfld (concld_fldn, 'fraction', pver, rad_data_avgflag,&
         'radiation input: convective cloud fraction',phys_decomp)
    call addfld (rel_fldn,    'micron',   pver, rad_data_avgflag,&
         'radiation input: effective liquid drop radius',phys_decomp)
    call addfld (rei_fldn,    'micron',   pver, rad_data_avgflag,&
         'radiation input: effective ice partical radius',phys_decomp)
    
    if (mg_microphys) then
       call addfld (dei_fldn,    'micron',   pver, rad_data_avgflag,&
            'radiation input: effective ice partical diameter',phys_decomp)
       call addfld (des_fldn,    'micron',   pver, rad_data_avgflag,&
            'radiation input: effective snow partical diameter',phys_decomp)
       call addfld (mu_fldn,     ' ',        pver, rad_data_avgflag,&
            'radiation input: ice gamma parameter for optics (radiation)',phys_decomp)
       call addfld (lambdac_fldn,' ',        pver, rad_data_avgflag,&
            'radiation input: slope of droplet distribution for optics (radiation)',phys_decomp)
       call addfld (iciwp_fldn,  'kg/m2',    pver, rad_data_avgflag,&
            'radiation input: In-cloud ice water path',phys_decomp)
       call addfld (iclwp_fldn,  'kg/m2',    pver, rad_data_avgflag,&
            'radiation input: In-cloud liquid water path',phys_decomp)
       call addfld (icswp_fldn,  'kg/m2',    pver, rad_data_avgflag,&
            'radiation input: In-cloud snow water path',phys_decomp)
       call addfld (rel_fn_fldn, 'microns',  pver, rad_data_avgflag,&
            'radiation input: ice effective drop size at fixed number (indirect effect)',phys_decomp)
       call addfld (cldfsnow_fldn, 'fraction', pver, rad_data_avgflag,&
            'radiation input: cloud liquid drops + snow',phys_decomp)
    endif

    call add_default (lndfrc_fldn,    rad_data_histfile_num, ' ')
    call add_default (icefrc_fldn,    rad_data_histfile_num, ' ')
    call add_default (snowh_fldn,     rad_data_histfile_num, ' ')
    call add_default (landm_fldn,     rad_data_histfile_num, ' ')
    call add_default (asdir_fldn,     rad_data_histfile_num, ' ')
    call add_default (asdif_fldn,     rad_data_histfile_num, ' ')
    call add_default (aldir_fldn,     rad_data_histfile_num, ' ')
    call add_default (aldif_fldn,     rad_data_histfile_num, ' ')

    call add_default (coszen_fldn,    rad_data_histfile_num, ' ')
    call add_default (asdir_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (asdif_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (aldir_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (aldif_pos_fldn, rad_data_histfile_num, ' ')

    call add_default (lwup_fldn,      rad_data_histfile_num, ' ')
    call add_default (ts_fldn,        rad_data_histfile_num, ' ')
    call add_default (temp_fldn,      rad_data_histfile_num, ' ')
    call add_default (pdel_fldn,      rad_data_histfile_num, ' ')
    call add_default (pdeldry_fldn,   rad_data_histfile_num, ' ')
    call add_default (pmid_fldn,      rad_data_histfile_num, ' ')
    call add_default (watice_fldn,    rad_data_histfile_num, ' ')
    call add_default (watliq_fldn,    rad_data_histfile_num, ' ')
    call add_default (watvap_fldn,    rad_data_histfile_num, ' ')
    call add_default (zint_fldn,      rad_data_histfile_num, ' ')
    call add_default (pint_fldn,      rad_data_histfile_num, ' ')

    call add_default (cld_fldn,       rad_data_histfile_num, ' ')
    call add_default (concld_fldn,    rad_data_histfile_num, ' ')
    call add_default (rel_fldn,       rad_data_histfile_num, ' ')
    call add_default (rei_fldn,       rad_data_histfile_num, ' ')
    
    if (mg_microphys) then
       call add_default (dei_fldn,       rad_data_histfile_num, ' ')
       call add_default (des_fldn,       rad_data_histfile_num, ' ')
       call add_default (mu_fldn,        rad_data_histfile_num, ' ')
       call add_default (lambdac_fldn,   rad_data_histfile_num, ' ')
       call add_default (iciwp_fldn,     rad_data_histfile_num, ' ')
       call add_default (iclwp_fldn,     rad_data_histfile_num, ' ')
       call add_default (icswp_fldn,     rad_data_histfile_num, ' ')
       call add_default (rel_fn_fldn,    rad_data_histfile_num, ' ')
       call add_default (cldfsnow_fldn,  rad_data_histfile_num, ' ')
    endif

    ! rad constituents

    call rad_cnst_get_info( 0, ngas=ngas, naero=naer )
    nrad_cnsts = ngas+naer
    allocate( rad_cnstnames(nrad_cnsts) )
    allocate( rad_cnstsources(nrad_cnsts) )
    allocate( rad_cnstindices(nrad_cnsts) )
    call rad_cnst_get_info( 0, gasnames=rad_cnstnames(1:ngas), aernames=rad_cnstnames(ngas+1:ngas+naer), &
                                 gassources=rad_cnstsources(1:ngas), aersources=rad_cnstsources(ngas+1:ngas+naer), &
                                 gasindices=rad_cnstindices(1:ngas), aerindices=rad_cnstindices(ngas+1:ngas+naer) )

    long_name_description = ' used in rad climate calculation'

    do i = 1, nrad_cnsts
       long_name = trim(rad_cnstnames(i))//' mass mixing ratio'//trim(long_name_description)
       name = 'rad_'//rad_cnstnames(i)
       call addfld(trim(name), 'kg/kg', pver, rad_data_avgflag, trim(long_name), phys_decomp)
       call add_default (trim(name), rad_data_histfile_num, ' ')
    end do

  end subroutine init_rad_data

  !================================================================================================
  !================================================================================================
  subroutine output_rad_data( pbuf, state, cam_in, landm, coszen )

    use physics_types,    only: physics_state
    use camsrfexch_types, only: cam_in_t     
    use phys_buffer,      only: pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx, pbuf_size_max
    use constituents,     only: cnst_get_ind

    implicit none

    type(pbuf_fld),      intent(in) :: pbuf(pbuf_size_max)
    type(physics_state), intent(in) :: state
    type(cam_in_t),      intent(in) :: cam_in
    real(r8),            intent(in) :: landm(pcols)
    real(r8),            intent(in) :: coszen(pcols)

    ! Local variables
    integer :: i
    character(len=1)  :: source
    character(len=32) :: name
    real(r8) :: mmr(pcols,pver)

    integer :: lchnk, itim, ifld
    integer :: ixcldice              ! cloud ice water index
    integer :: ixcldliq              ! cloud liquid water index
    integer :: icol
    integer :: ncol

    ! surface albedoes weighted by (positive cosine zenith angle)
    real(r8):: coszrs_pos(pcols)    ! = max(coszrs,0)
    real(r8):: asdir_pos (pcols)    !
    real(r8):: asdif_pos (pcols)    !
    real(r8):: aldir_pos (pcols)    !
    real(r8):: aldif_pos (pcols)    !

    real(r8), pointer, dimension(:,:)  :: ptr

    if (.not.rad_data_output) return

    ! get index of (liquid+ice) cloud water
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

    lchnk = state%lchnk
    ncol = state%ncol

    do icol = 1, ncol
       coszrs_pos(icol)  = max(coszen(icol),0._r8)
    enddo
    asdir_pos(:ncol)  = cam_in%asdir(:ncol) * coszrs_pos(:ncol)
    asdif_pos(:ncol)  = cam_in%asdif(:ncol) * coszrs_pos(:ncol)
    aldir_pos(:ncol)  = cam_in%aldir(:ncol) * coszrs_pos(:ncol)
    aldif_pos(:ncol)  = cam_in%aldif(:ncol) * coszrs_pos(:ncol)

    call outfld(lndfrc_fldn, cam_in%landfrac,  pcols, lchnk)
    call outfld(icefrc_fldn, cam_in%icefrac,   pcols, lchnk)
    call outfld(snowh_fldn,  cam_in%snowhland, pcols, lchnk)
    call outfld(landm_fldn,  landm,            pcols, lchnk)
    call outfld(temp_fldn,   state%t,               pcols, lchnk   )
    call outfld(pdel_fldn,   state%pdel,            pcols, lchnk   )
    call outfld(pdeldry_fldn,state%pdeldry,         pcols, lchnk   )
    call outfld(watice_fldn, state%q(:,:,ixcldice), pcols, lchnk   )
    call outfld(watliq_fldn, state%q(:,:,ixcldliq), pcols, lchnk   )
    call outfld(watvap_fldn, state%q(:,:,1),        pcols, lchnk   )
    call outfld(zint_fldn,   state%zi,              pcols, lchnk   )
    call outfld(pint_fldn,   state%pint,            pcols, lchnk   )
    call outfld(pmid_fldn,   state%pmid,            pcols, lchnk   )

    call outfld(asdir_fldn, cam_in%asdir, pcols, lchnk   )
    call outfld(asdif_fldn, cam_in%asdif, pcols, lchnk   )
    call outfld(aldir_fldn, cam_in%aldir, pcols, lchnk   )
    call outfld(aldif_fldn, cam_in%aldif, pcols, lchnk   )

    call outfld(coszen_fldn, coszrs_pos, pcols, lchnk   )
    call outfld(asdir_pos_fldn, asdir_pos, pcols, lchnk   )
    call outfld(asdif_pos_fldn, asdif_pos, pcols, lchnk   )
    call outfld(aldir_pos_fldn, aldir_pos, pcols, lchnk   )
    call outfld(aldif_pos_fldn, aldif_pos, pcols, lchnk   )

    call outfld(lwup_fldn,  cam_in%lwup,  pcols, lchnk   )
    call outfld(ts_fldn,    cam_in%ts,    pcols, lchnk   )

    itim = pbuf_old_tim_idx()

    ptr => pbuf(cld_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
    call outfld(cld_fldn,    ptr,    pcols, lchnk   )
    ptr => pbuf(concld_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
    call outfld(concld_fldn, ptr, pcols, lchnk   )
    ptr  => pbuf(rel_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
    call outfld(rel_fldn,    ptr,    pcols, lchnk   )
    ptr  => pbuf(rei_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
    call outfld(rei_fldn,    ptr,    pcols, lchnk   )

    if (mg_microphys) then
       ptr => pbuf(dei_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
       call outfld(dei_fldn,    ptr,    pcols, lchnk   )       
       ptr => pbuf(des_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
       call outfld(des_fldn,    ptr,    pcols, lchnk   )       
       ptr => pbuf(mu_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
       call outfld(mu_fldn,    ptr,    pcols, lchnk   ) 
       ptr => pbuf(lambdac_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
       call outfld(lambdac_fldn,    ptr,    pcols, lchnk   )       
       ptr => pbuf(iciwp_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
       call outfld(iciwp_fldn,    ptr,    pcols, lchnk   )       
       ptr => pbuf(iclwp_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
       call outfld(iclwp_fldn,    ptr,    pcols, lchnk   )       
       ptr => pbuf(icswp_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
       call outfld(icswp_fldn,    ptr,    pcols, lchnk   )       
       ptr => pbuf(rel_fn_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
       call outfld(rel_fn_fldn,    ptr,    pcols, lchnk   )       
       ptr => pbuf(cldfsnow_ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
       call outfld(cldfsnow_fldn,    ptr,    pcols, lchnk   )
    endif

    ! output mixing ratio of rad constituents 

    do i = 1, nrad_cnsts

       name = 'rad_'//rad_cnstnames(i)
       source = rad_cnstsources(i)

       select case( source )
       case ('P')
          mmr(:ncol,:) = state%q(:ncol,:,rad_cnstindices(i))
       case ('D')
          mmr(:ncol,:) = pbuf(rad_cnstindices(i))%fld_ptr(1,:ncol,:,lchnk,1)
       end select
       call outfld(trim(name), mmr(:ncol,:), ncol, lchnk)

    end do

  end subroutine output_rad_data


end module radiation_data
