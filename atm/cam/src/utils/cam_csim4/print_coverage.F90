
subroutine print_coverage (string1, string2, variable, factor)

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,    only: pcols, begchunk, endchunk
   use pmgrid,    only: plon, plat
   use ice_spmd  ,only: masterproc
   use commap,    only: w
   use phys_grid, only: gather_chunk_to_field, get_ncols_p
   use rgrid,     only: nlon
   use physconst, only: pi,rearth

   implicit none
!
! Arguments
!
   character(len=*), intent(in) :: string1
   character(len=*), intent(in) :: string2
   real(r8), intent(in) :: variable(pcols,begchunk:endchunk)
   real(r8), intent(in) :: factor
!
! Local workspace
!
   real(r8) :: hemis_area ! = 2._r8*pi*rearth*rearth
   real(r8) :: var_field(plon,plat)
   real(r8) :: shcoverage, nhcoverage
   real(r8) :: totsum
   real(r8) :: sum

   integer :: i,j

   hemis_area = 2._r8*pi*rearth*rearth

   call gather_chunk_to_field (1, 1, 1, plon, variable, var_field)
   
   if (masterproc) then
      totsum = 0._r8
      do j=1,plat/2
         sum = 0._r8
         do i=1,nlon(j)
            sum = sum + var_field(i,j)*w(j)
         end do
         totsum = totsum + sum/nlon(j)
      end do
!
! Assume the weights w sum to 1 for the hemisphere
!
      shcoverage = totsum*hemis_area
      shcoverage = shcoverage*factor        ! conversion factor
      write(6,'(a,a,f9.4,a)')string1,' shcoverage=', shcoverage, string2

      totsum = 0._r8
      do j=plat/2+1,plat
         sum = 0._r8
         do i=1,nlon(j)
            sum = sum + var_field(i,j)*w(j)
         end do
         totsum = totsum + sum/nlon(j)
      end do
!
! Assume the weights w sum to 1 for the hemisphere
!
      nhcoverage = totsum*hemis_area
      nhcoverage = nhcoverage*factor        ! conversion factor
      write(6,'(a,a,f9.4,a)')string1,' nhcoverage=', nhcoverage, string2

   end if
   return
end subroutine print_coverage
