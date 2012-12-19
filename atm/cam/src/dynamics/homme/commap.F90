module commap

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plon, plev, plat

   real(r8), pointer :: w(:)            ! gaussian weights (hemisphere)

   real(r8), pointer:: clat(:)         ! model latitudes (radians)
   real(r8), pointer:: clon(:,:)   ! model longitudes (radians)
   real(r8), pointer:: latdeg(:)       ! model latitudes (degrees)

   real(r8), pointer:: londeg(:,:) ! model longitudes (degrees)

end module commap
