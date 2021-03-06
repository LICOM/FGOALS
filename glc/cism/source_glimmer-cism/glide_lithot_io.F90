!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WARNING: this file was automatically generated on
! Tue, 15 Jun 2010 22:12:39 +0000
! from ncdf_template.F90.in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  ncdf_template.f90 - part of the Glimmer_CISM ice model   + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010
! Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
! This file is part of Glimmer-CISM.
!
! Glimmer-CISM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or (at
! your option) any later version.
!
! Glimmer-CISM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
!
! Glimmer-CISM is hosted on BerliOS.de:
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define NCO outfile%nc
#define NCI infile%nc


module glide_lithot_io
  !*FD template for creating subsystem specific I/O routines
  !*FD written by Magnus Hagdorn, 2004

  private :: get_xtype

  character(len=*),private,parameter :: hotvars = ' litho_temp '

contains

  !*****************************************************************************
  ! netCDF output
  !*****************************************************************************
  subroutine glide_lithot_io_createall(model,data,outfiles)
    !*FD open all netCDF files for output
    use glide_types
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none
    type(glide_global_type) :: model
    type(glide_global_type), optional :: data
    type(glimmer_nc_output),optional,pointer :: outfiles
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    do while(associated(oc))
       if (present(data)) then
          call glide_lithot_io_create(oc,model,data)
       else
          call glide_lithot_io_create(oc,model)
       end if
       oc=>oc%next
    end do
  end subroutine glide_lithot_io_createall

  subroutine glide_lithot_io_writeall(data,model,atend,outfiles,time)
    !*FD if necessary write to netCDF files
    use glide_types
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none
    type(glide_global_type) :: data
    type(glide_global_type) :: model
    logical, optional :: atend
    type(glimmer_nc_output),optional,pointer :: outfiles
    real(sp),optional :: time

    ! local variables
    type(glimmer_nc_output), pointer :: oc
    logical :: forcewrite=.false.

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    if (present(atend)) then
       forcewrite = atend
    end if

    do while(associated(oc))
#ifdef HAVE_AVG
       if (oc%do_averages) then
          call glide_lithot_avg_accumulate(oc,data,model)
       end if
#endif
       call glimmer_nc_checkwrite(oc,model,forcewrite,time)
       if (oc%nc%just_processed) then
          ! write standard variables
          call glide_lithot_io_write(oc,data)
#ifdef HAVE_AVG
          if (oc%do_averages) then
             call glide_lithot_avg_reset(oc,data)
          end if
#endif
       end if
       oc=>oc%next
    end do
  end subroutine glide_lithot_io_writeall
  
  subroutine glide_lithot_io_create(outfile,model,data)
    use glide_types
    use glide_types
    use glimmer_ncdf
    use glimmer_map_types
    use glimmer_log
    use glimmer_paramets
    use glimmer_scales
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    type(glide_global_type) :: model
    type(glide_global_type), optional :: data

    integer status,varid,pos

    integer :: lithoz_dimid
    integer :: time_dimid
    integer :: x1_dimid
    integer :: y1_dimid

    ! defining dimensions
    if (.not.outfile%append) then
       status = nf90_def_dim(NCO%id,'lithoz',model%lithot%nlayer,lithoz_dimid)
    else
       status = nf90_inq_dimid(NCO%id,'lithoz',lithoz_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inq_dimid(NCO%id,'time',time_dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inq_dimid(NCO%id,'x1',x1_dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inq_dimid(NCO%id,'y1',y1_dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)

    NCO%vars = ' '//trim(NCO%vars)//' '
    ! expanding hotstart variables
    pos = index(NCO%vars,' hot ') 
    if (pos.ne.0) then
       NCO%vars = NCO%vars(:pos)//NCO%vars(pos+4:)
       NCO%hotstart = .true.
    end if
    if (NCO%hotstart) then
       NCO%vars = trim(NCO%vars)//hotvars
    end if
    ! checking if we need to handle time averages
    pos = index(NCO%vars,"_tavg")
    if (pos.ne.0) then
       outfile%do_averages = .True.
    end if    

    !     lithoz -- vertical coordinate of lithosphere layer
    if (.not.outfile%append) then
       call write_log('Creating variable lithoz')
       status = nf90_def_var(NCO%id,'lithoz',get_xtype(outfile,NF90_FLOAT),(/lithoz_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NCO%id, varid, 'long_name', 'vertical coordinate of lithosphere layer')
       status = nf90_put_att(NCO%id, varid, 'units', 'meter')
     end if

    !     litho_temp -- lithosphere temperature
    pos = index(NCO%vars,' litho_temp ')
    status = nf90_inq_varid(NCO%id,'litho_temp',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
       call write_log('Creating variable litho_temp')
       status = nf90_def_var(NCO%id,'litho_temp',get_xtype(outfile,NF90_FLOAT),(/x1_dimid, y1_dimid, lithoz_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NCO%id, varid, 'long_name', 'lithosphere temperature')
       status = nf90_put_att(NCO%id, varid, 'units', 'degree_Celsius')
       if (glimmap_allocated(model%projection)) then
          status = nf90_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = nf90_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     end if

  end subroutine glide_lithot_io_create

  subroutine glide_lithot_io_write(outfile,data)
    use glide_types
    use glimmer_ncdf
    use glimmer_paramets
    use glimmer_scales
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: data
    !*FD the model instance

    ! local variables
    real tavgf
    integer status, varid
    integer up
     
    tavgf = outfile%total_time
    if (tavgf.ne.0.) then
       tavgf = 1./tavgf
    end if

    ! write variables
    status = nf90_inq_varid(NCO%id,'litho_temp',varid)
    if (status .eq. nf90_noerr) then
       status = nf90_put_var(NCO%id, varid, &
            data%lithot%temp, (/1,1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

  end subroutine glide_lithot_io_write

  !*****************************************************************************
  ! netCDF input
  !*****************************************************************************  
  subroutine glide_lithot_io_readall(data,model)
    !*FD read from netCDF file
    use glide_types
    use glide_types
    use glimmer_ncio
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: data
    type(glide_global_type) :: model

    ! local variables
    type(glimmer_nc_input), pointer :: ic    

    ic=>model%funits%in_first
    do while(associated(ic))
       call glimmer_nc_checkread(ic,model)
       if (ic%nc%just_processed) then
          call glide_lithot_io_read(ic,data)
       end if
       ic=>ic%next
    end do
  end subroutine glide_lithot_io_readall

  subroutine glide_lithot_io_read(infile,data)
    !*FD read variables from a netCDF file
    use glimmer_log
    use glimmer_ncdf
    use glide_types
    use glimmer_paramets
    use glimmer_scales
    implicit none
    type(glimmer_nc_input), pointer :: infile
    !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: data
    !*FD the model instance

    ! local variables
    integer status,varid
    integer up
    real(dp) :: scaling_factor

    ! read variables
    status = nf90_inq_varid(NCI%id,'litho_temp',varid)
    if (status .eq. nf90_noerr) then
       call write_log('  Loading litho_temp')
       status = nf90_get_var(NCI%id, varid, &
            data%lithot%temp, (/1,1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling litho_temp",GM_DIAGNOSTIC)
          data%lithot%temp = data%lithot%temp*scaling_factor
       end if
    end if

  end subroutine glide_lithot_io_read

  subroutine glide_lithot_io_checkdim(infile,model,data)
    !*FD check if dimension sizes in file match dims of model
    use glimmer_log
    use glimmer_ncdf
    use glide_types
    use glide_types
    implicit none
    type(glimmer_nc_input), pointer :: infile
    !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: model
    type(glide_global_type), optional :: data

    integer status,dimid,dimsize
    character(len=150) message

    ! check dimensions
    status = nf90_inq_dimid(NCI%id,'lithoz',dimid)
    if (dimid.gt.0) then
       status = nf90_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%lithot%nlayer) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size lithoz does not match: ', &
               model%lithot%nlayer
          call write_log(message,GM_FATAL)
       end if
    end if
  end subroutine glide_lithot_io_checkdim

  !*****************************************************************************
  ! calculating time averages
  !*****************************************************************************  
#ifdef HAVE_AVG
  subroutine glide_lithot_avg_accumulate(outfile,data,model)
    use glide_types
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: model
    type(glide_global_type) :: data

    ! local variables
    real :: factor
    integer status, varid

    ! increase total time
    outfile%total_time = outfile%total_time + model%numerics%tinc
    factor = model%numerics%tinc

  end subroutine glide_lithot_avg_accumulate

  subroutine glide_lithot_avg_reset(outfile,data)
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: data

    ! local variables
    integer status, varid

    ! reset total time
    outfile%total_time = 0.

  end subroutine glide_lithot_avg_reset
#endif

  !*********************************************************************
  ! some private procedures
  !*********************************************************************

  !> apply default type to be used in netCDF file
  integer function get_xtype(outfile,xtype)
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_output), pointer :: outfile !< derived type holding information about output file
    integer, intent(in) :: xtype                !< the external netCDF type

    get_xtype = xtype
    
    if (xtype.eq.NF90_REAL .and. outfile%default_xtype.eq.NF90_DOUBLE) then
       get_xtype = NF90_DOUBLE
    end if
    if (xtype.eq.NF90_DOUBLE .and. outfile%default_xtype.eq.NF90_REAL) then
       get_xtype = NF90_REAL
    end if
  end function get_xtype

  !*********************************************************************
  ! lots of accessor subroutines follow
  !*********************************************************************

end module glide_lithot_io
