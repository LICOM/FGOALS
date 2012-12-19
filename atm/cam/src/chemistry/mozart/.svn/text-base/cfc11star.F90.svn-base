!---------------------------------------------------------------------------------
! Manages the CFC11* for radiation 
!  4 Dec 2009 -- Francis Vitt created
!---------------------------------------------------------------------------------
module cfc11star

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  use phys_buffer,  only : pbuf, pbuf_add
  use abortutils,   only : endrun
  use ppgrid,       only : pver, begchunk, endchunk
  use spmd_utils,   only : masterproc

  implicit none
  save 

  private
  public :: register_cfc11star
  public :: update_cfc11star
  public :: init_cfc11star

  logical :: do_cfc11star
  character(len=16), parameter :: pbufname = 'CFC11STAR'
  integer :: pbf_idx
 
  integer, pointer :: cfc11_ndx
  integer, pointer :: cfc113_ndx
  integer, pointer :: ccl4_ndx
  integer, pointer :: ch3ccl3_ndx
  integer, pointer :: hcfc22_ndx
  integer, pointer :: cf2clbr_ndx
  integer, pointer :: cf3br_ndx
  integer, target :: indices(7)
  
  real(r8) :: rel_rf(7)

contains

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
  subroutine register_cfc11star

    implicit none
    
    call pbuf_add(pbufname,'global',1,pver,1,pbf_idx)

  endsubroutine register_cfc11star

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
  subroutine init_cfc11star
    use constituents, only : cnst_get_ind
    use cam_history,  only : addfld, phys_decomp
    use infnan,       only : nan

    implicit none

    real(r8), parameter :: cfc_rf(7)  =  (/ 0.25, 0.30, 0.13, 0.06, 0.20, 0.30, 0.32 /)    ! W/m2/ppb

    cfc11_ndx   => indices(1)
    cfc113_ndx  => indices(2)
    ccl4_ndx    => indices(3)
    ch3ccl3_ndx => indices(4)
    hcfc22_ndx  => indices(5)
    cf2clbr_ndx => indices(6)
    cf3br_ndx   => indices(7)
    
    call cnst_get_ind('CFC11',  cfc11_ndx,   abort=.false.)
    call cnst_get_ind('CFC113', cfc113_ndx,  abort=.false.)
    call cnst_get_ind('CCL4',   ccl4_ndx,    abort=.false.)
    call cnst_get_ind('CH3CCL3',ch3ccl3_ndx, abort=.false.)
    call cnst_get_ind('HCFC22', hcfc22_ndx,  abort=.false.)
    call cnst_get_ind('CF2CLBR',cf2clbr_ndx, abort=.false.)
    call cnst_get_ind('CF3BR',  cf3br_ndx,   abort=.false.)

    do_cfc11star = all(indices(:)>0)

    if (.not.do_cfc11star) return

    pbuf(pbf_idx)%fld_ptr = nan

    rel_rf(:) = cfc_rf(:) / cfc_rf(1)
    call addfld(pbufname,'kg/kg',pver,'A','cfc11star for radiation', phys_decomp )
    
    if (masterproc) then
       write(iulog,*) 'init_cfc11star: CFC11STAR is added to pbuf for radiation'
    endif
  end subroutine init_cfc11star

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
  subroutine update_cfc11star( phys_state )
    use cam_history,  only : outfld
    use physics_types,only : physics_state

    implicit none

    type(physics_state), intent(in):: phys_state(begchunk:endchunk)                 

    integer :: lchnk, ncol
    integer :: c
    real(r8), pointer :: cf11star(:,:)

    if (.not.do_cfc11star) return
    
    do c = begchunk,endchunk
       lchnk = phys_state(c)%lchnk
       ncol = phys_state(c)%ncol

       cf11star => pbuf(pbf_idx)%fld_ptr(1,:,:,lchnk,1)

       cf11star(:ncol,:) = &
            phys_state(c)%q(:ncol,:,cfc11_ndx)  * rel_rf(1) + &
            phys_state(c)%q(:ncol,:,cfc113_ndx) * rel_rf(2) + &
            phys_state(c)%q(:ncol,:,ccl4_ndx)   * rel_rf(3) + &
            phys_state(c)%q(:ncol,:,ch3ccl3_ndx)* rel_rf(4) + &
            phys_state(c)%q(:ncol,:,hcfc22_ndx) * rel_rf(5) + &
            phys_state(c)%q(:ncol,:,cf2clbr_ndx)* rel_rf(6) + &
            phys_state(c)%q(:ncol,:,cf3br_ndx)  * rel_rf(7)

       call outfld( pbufname, cf11star(:ncol,:), ncol, lchnk )

    enddo

  endsubroutine update_cfc11star

endmodule cfc11star
