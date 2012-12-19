!---------------------------------------------------------------------
! Manages the storage of non-transported short-lived chemical species
! in the physics buffer.
!
! Created by: Francis Vitt -- 20 Aug 2008
!---------------------------------------------------------------------
module short_lived_species

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods,    only : slvd_lst, nslvd, gas_pcnst
  use cam_logfile,  only : iulog
  use ppgrid,       only : pcols, pver, begchunk, endchunk
  use spmd_utils,   only : masterproc
  use phys_buffer,  only : pbuf

  implicit none

  save
  private
  public :: map
  public :: register_short_lived_species
  public :: initialize_short_lived_species
  public :: set_short_lived_species
  public :: get_short_lived_species
  public :: slvd_index
  public :: pbf_idx

  integer :: pbf_idx
  integer :: map(nslvd)

  character(len=16), parameter :: pbufname = 'ShortLivedSpecies'

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine register_short_lived_species
    use phys_buffer, only : pbuf_add

    implicit none

    integer :: m

    if ( nslvd < 1 ) return

    call pbuf_add(pbufname,'global',1,pver,nslvd,pbf_idx)

  end subroutine register_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine initialize_short_lived_species(ncid_ini)
    use ioFileMod,      only : getfil
    use error_messages, only : handle_ncerr
    use dycore,         only : dycore_is
    use mo_tracname,    only : solsym
    use ncdio_atm,      only : infld
    use infnan,         only : nan
    use pio,            only : file_desc_t
    implicit none

    type(file_desc_t), intent(inout) :: ncid_ini

    integer :: m,n
    character(len=8) :: fieldname
    character(len=4) :: dim1name
    logical :: found

    if ( nslvd < 1 ) return

    found = .false.

    if(dycore_is('homme')) then  
       dim1name='ncol'
    else
       dim1name='lon'
    end if

    pbuf(pbf_idx)%fld_ptr(:,:,:,:,:) = nan

    do m=1,nslvd
       n = map(m)
       fieldname = solsym(n)
       call infld( fieldname,ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                   pbuf(pbf_idx)%fld_ptr(1,:,:,:,m), found, grid_map='PHYS')
       if(.not. found) then
          pbuf(pbf_idx)%fld_ptr(:,:,:,:,m) = 0._r8
          if (masterproc) write(iulog,*) 'short_lived_species: '//trim(fieldname)//' initialized to 0.'
       end if
    enddo

  endsubroutine initialize_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine set_short_lived_species( q, lchnk, ncol )
    implicit none 

    real(r8), intent(in) :: q(pcols,pver,gas_pcnst)
    integer,  intent(in) :: lchnk, ncol

    integer :: m,n

    if ( nslvd < 1 ) return

    do m=1,nslvd
       n = map(m)
       pbuf(pbf_idx)%fld_ptr(1,:ncol,:,lchnk,m) = q(:ncol,:,n)
    enddo

  endsubroutine set_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine get_short_lived_species( q, lchnk, ncol )
    implicit none 

    real(r8), intent(inout) :: q(pcols,pver,gas_pcnst)
    integer,  intent(in) :: lchnk, ncol

    integer :: m,n 

    if ( nslvd < 1 ) return

    do m=1,nslvd
       n = map(m)
       q(:ncol,:,n) = pbuf(pbf_idx)%fld_ptr(1,:ncol,:,lchnk,m)
    enddo

  endsubroutine get_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  function slvd_index( name )
    implicit none

    character(len=*) :: name
    integer :: slvd_index

    integer :: m

    slvd_index = -1

    if ( nslvd < 1 ) return

    do m=1,nslvd
       if ( name == slvd_lst(m) ) then
          slvd_index = m
          return 
       endif
    enddo

  endfunction slvd_index

end module short_lived_species
