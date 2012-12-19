!
! test dshr_tInterp as well as dshr_bundle and dshr_domain code
!
program test_dshr_tInterp
   use shr_sys_mod
   use shr_kind_mod
   use dshr_kind
   use test_mod
   use dshr_tInterp
   use bundle_expected
   use shr_const_mod
   use dshr_bundle
   use dshr_domain
   use dshr_tInterp
   use shr_string_mod
   use shr_orb_mod
   use shr_date_mod

   implicit none

   real(R8) :: f1, f2
   character(SHR_KIND_CS) :: alogo, filltype, dom_name, name
   character(SHR_KIND_CX) :: fieldList, fieldList2, fieldListBound
   real(r8), parameter :: dtime = 3600
   integer :: dt
   integer, parameter :: LIN_TEST = 1, LOWER_TEST = 2, UPPER_TEST = 3, &
                         NEAREST_TEST = 4, COSZ_TEST=5, num_tests = 5, num_times = 6*3600/int(dtime), &
                         num_types = 4, num_fail = 7, num_periods = 4
   integer, parameter :: INDX_FILL=1, SINCOS_FILL=2, NFIELD_FILL=3, COSZ_FILL=4
   real(R8) :: factors(2), mfactor(num_periods+1)
   integer :: ntimes
   integer :: date_lb(num_periods+1), date_ub(num_periods+2), date_in
   integer :: sec_lb(num_periods+1), sec_ub(num_periods+2), sec_in, sec_ub_eff
   type(dshr_bundle_bundleType) :: bunlb         ! lower bound
   type(dshr_bundle_bundleType) :: bunub         ! upper bound
   type(dshr_bundle_bundleType) :: bun           ! output bundle
   type(dshr_bundle_bundleType) :: expected_bun  ! expected bundle
   type(dshr_domain_domainType), target :: domain, domain2
   type(dshr_domain_domainType), pointer :: domainP1, domainP2, domainP3
   type(shr_date) :: sdate, sdateub, sdatenext
   integer :: n, i, ni, nj, nf, rc, type, date_start, sec_start, b, t, kfld, kfldLB
   logical :: status
   real(r8), pointer :: data(:,:,:), data2(:,:,:), lat(:,:), lon(:,:), data3(:,:,:)
   real(r8) :: scale
   real(r8) :: orb_eccen, orb_obliq, orb_mvelp,     &
               orb_obliqr, orb_lambm0 , orb_mvelpp, declin, eccf, calday
   integer :: orb_iyear_AD

   call test_init( num_types*num_tests*(num_periods+1)*(num_times+1) + num_fail + 1 - num_types*num_tests )
   !
   ! Create a domain
   !
   dom_name = "test_domain"
   ni = 40
   nj = 75
   call dshr_domain_fill(domain,dom_name,ni,nj)
   allocate( lat(ni,nj) )
   allocate( lon(ni,nj) )
   !
   ! Loop over the different types of fill type to fill the bundles up with
   !
   do type = 1, num_types

      if (      type == INDX_FILL )then
         filltype = "index"
         fieldList = "T:U:V:W:TS:VAL:FSDS" 
      else if ( type == SINCOS_FILL )then
         filltype = "sincos"
         fieldList = "T:FSDS:U"
      else if ( type == NFIELD_FILL )then
         filltype = "nfield"
         fieldList = "T:U:V:TS:PRECT:FSDS:FLDS"
      else if ( type == COSZ_FILL )then
         filltype = "cosz"
         fieldList = "FSDS:T:U:V"
         scale = 1000.0_r8
      else
         call shr_sys_abort( "type is out of bounds" )
      end if

      kfld = shr_string_listGetIndexF( fieldList, "FSDS" )

      write(*,*) "fill type = ", trim(filltype)
      
      !
      ! Now test the actual interpolations
      !
      do n = 1, num_tests

         fieldListBound = trim(fieldList)
         if (      n == LIN_TEST     )then
            alogo = 'linear'
         else if ( n == LOWER_TEST   )then
            alogo = 'lower'
         else if ( n == UPPER_TEST   )then
            alogo = 'upper'
         else if ( n == NEAREST_TEST )then
            alogo = 'nearest'
         else if ( n == COSZ_TEST )then
            fieldListBound = "FSDS"//trim(dshr_tInterp_cosZScaleList)
            alogo = 'coszen'
            orb_iyear_AD = 1990
            call shr_orb_params( orb_iyear_AD, orb_eccen, orb_obliq, orb_mvelp,     &
           &                     orb_obliqr, orb_lambm0 , orb_mvelpp, log_print=.true. )
         else
            alogo = 'undefined'
         end if
         nf = shr_string_listGetNum( fieldList )
         ! Create and fill the bundles
         name = filltype
         call dshr_bundle_create( bunlb,        domain, "lb_"//name, fieldListBound )
         call dshr_bundle_create( bunub,        domain, "ub_"//name, fieldListBound )
         call dshr_bundle_create( bun,          domain, "bun_"//name, fieldList )
         call dshr_bundle_create( expected_bun, domain, "expected_"//name, fieldList )
         kfldLB = shr_string_listGetIndexF( fieldListBound, "FSDS" )
      
         write(*,*) "Test type: ", trim(alogo)

         ! 4 periods per day (total to one day) + 1 period into the next day
         mfactor = (/ 2.0_r8, 6.0_r8, 4.0_r8, 1.0_r8, 1.5_r8 /)
         date_lb = 19481231
         sec_lb  = (/ 0, 3600*6, 3600*12, 3600*18, &
                      nint(0.999_r8 * 3600._r8 * 24.0_r8) /)
         date_ub(1:num_periods) = date_lb(1)
         date_ub(num_periods+1:)= 19490101
         sec_ub (1:num_periods) = sec_lb (2:num_periods+1)
         sec_ub (num_periods+1) = 3600*6
         sec_ub (num_periods+2) = 3600*12
         do b = 1, num_periods+1
            sdateub = shr_date_initCDate( date_ub(b), num_periods, sec_ub(b) )
            sdatenext = shr_date_initCDate( date_ub(b+1), num_periods, sec_ub(b+1) )
            call dshr_bundle_putDate( bunlb, cdate=date_lb(b), sec=sec_lb(b) )
            call dshr_bundle_putDate( bunub, cdate=date_ub(b), sec=sec_ub(b) )
            if ( date_ub(b) > date_lb(b) )then
               sec_ub_eff = sec_ub(b) + 86400
            else
               sec_ub_eff = sec_ub(b)
            end if

            ! Fill bundles for lower and upper bounds
            if ( b == 1 )then
               if ( type == COSZ_FILL )then
                  call dshr_bundle_fill( bunlb, type="sincos" )
                  call dshr_bundle_fill( bunub, type="sincos" )
                  call dshr_bundle_domainPtr( bunlb, domainP1 )
                  call bundle_fill_cosz( scale, orb_eccen, orb_mvelpp, &
                                         orb_lambm0, orb_obliqr, sdateub, domainP1, &
                                         bunlb, kfld=kfldLB)
                  nullify( domainP1 )
                  call dshr_bundle_domainPtr( bunlb, domainP1 )
                  call bundle_fill_cosz( scale, orb_eccen, orb_mvelpp, &
                                         orb_lambm0, orb_obliqr, sdatenext, domainP1, &
                                         bunub, kfld=kfldLB )
                  nullify( domainP1 )
               else
                  call dshr_bundle_fill( bunlb, type=filltype )
                  call dshr_bundle_fill( bunub, type=filltype )
                  ! Multiply upper bound by a factor
                  call dshr_bundle_assignPtr( bunub, data )
                  data = data*mfactor(b)
                  nullify( data )
               end if
            else
               call dshr_bundle_assignPtr( bunub, data )
               call dshr_bundle_assignPtr( bunlb, data2 )
               data2 = data   ! Copy lower bound from upper bound data
               nullify(data)
               nullify(data2)
               if ( type /= COSZ_FILL )then
                  call dshr_bundle_fill( bunub, type=filltype )
                  ! Multiply upper bound by a factor
                  call dshr_bundle_assignPtr( bunub, data )
                  data = data*mfactor(b)
                  nullify( data )
               else
                  call dshr_bundle_fill( bunub, type="sincos")
                  call dshr_bundle_domainPtr( bunlb, domainP1 )
                  call bundle_fill_cosz( scale, orb_eccen, orb_mvelpp, &
                                         orb_lambm0, orb_obliqr, sdatenext, domainP1, &
                                         bunub, kfld=kfldLB )
                  nullify( domainP1 )
               end if
            end if
            if ( n == COSZ_TEST )then
               call dshr_tInterp_calcCosZenNorm( bunlb, sdateUB, orb_eccen,  &
                                                 orb_mvelpp, orb_lambm0, orb_obliqr )
            end if

            ! 6 hour blocks
            date_in = date_lb(b)
            sec_in  = sec_lb(b)
            sdate = shr_date_initCDate( date_in, num_periods*num_times, sec_in )
         
            do i = 1, num_times+1
               if ( sdate > sdateUB ) cycle
               call shr_date_getCDate( sdate, date_in, sec_in )
               write(*,*) "date in ", date_in, "seconds in ", sec_in
               call dshr_bundle_putDate( bun,          cdate=date_in, sec=sec_in )
               call dshr_bundle_putDate( expected_bun, cdate=date_in, sec=sec_in )
               if ( date_in > date_lb(b) ) sec_in = sec_in + 86400
               if (      n == LIN_TEST     )then
                  f1 = real((sec_in - sec_lb(b)),SHR_KIND_R8) / &
                       real((sec_ub_eff - sec_lb(b)),SHR_KIND_R8)
                  factors = (/ 1.0_SHR_KIND_R8 - f1, f1 /)
                  call bundle_linear_alg( bunlb, bunub, factors, expected_bun )
               else if ( n == COSZ_TEST    )then
                  call dshr_bundle_assignPtr( bunlb, data2 )
                  date_start = date_lb(b)
                  sec_start  = sec_lb(b)
                  call dshr_bundle_domainPtr( bunlb, domainP1 )
                  call dshr_domain_getData( domainP1, lat, "lat" )
                  call dshr_domain_getData( domainP1, lon, "lon" )
                  nullify( domainP1 )
                  lat = lat * SHR_CONST_PI / 180._r8
                  lon = lon * SHR_CONST_PI / 180._r8
                  calday = shr_date_getJulian( sdate )
                  call shr_orb_decl(calday, orb_eccen, orb_mvelpp, orb_lambm0, &
                                    orb_obliqr, declin, eccf)
                  call dshr_bundle_fill( bun,          type="sincos" )
                  if ( type == COSZ_FILL )then
                     call dshr_bundle_domainPtr( expected_bun, domainP1 )
                     call dshr_bundle_fill( expected_bun, type="sincos" )
                     call bundle_fill_cosz( scale, orb_eccen, orb_mvelpp, &
                                            orb_lambm0, orb_obliqr, sdate, domainP1, &
                                            expected_bun, kfld=kfld )
                     nullify( domainP1 )
                  else
                     call dshr_bundle_assignPtr( expected_bun,   data3  )
                     call dshr_bundle_fill( expected_bun, type="sincos" )
                     dt = int(dtime)/45
                     ntimes = (sec_ub_eff - sec_lb(b))/dt
                     call scalebycosz( ni, nj, nf, ntimes, dt, data2(:,:,kfldLB), &
                                       lat, lon, calday, declin, date_start, &
                                       sec_start, data3, kfld, option="simpsons" )
                     nullify( data3 )
                  end if
                  nullify( data2 )
               else if ( n == LOWER_TEST   )then
                  call dshr_bundle_assignPtr( expected_bun,   data  )
                  call dshr_bundle_assignPtr( bunlb,          data2 )
                  data = data2
                  nullify( data )
                  nullify( data2 )
               else if ( n == UPPER_TEST   )then
                  call dshr_bundle_assignPtr( expected_bun,   data  )
                  call dshr_bundle_assignPtr( bunub,          data2 )
                  data = data2
                  nullify( data )
                  nullify( data2 )
               else if ( n == NEAREST_TEST )then
                  if ( sec_in-sec_lb(b) <= (sec_ub_eff-sec_lb(b))/2 )then
                     call dshr_bundle_assignPtr( bunlb, data2 )
                  else
                     call dshr_bundle_assignPtr( bunub, data2 )
                  end if
                  call dshr_bundle_assignPtr( expected_bun,   data  )
                  data = data2
                  nullify( data )
                  nullify( data2 )
               end if
               call dshr_tInterp_data(bunlb,bunub,bun,algo=alogo)
               if (      n == LIN_TEST     )then
                  status = bundle_closeto_expected( bun, expected_bun, eps=1.e-4_r8 )
                  call test_is( status, "Test if resultant bundle is close to what is expected :: "//trim(name)//" "//trim(alogo) )
               else if ( n == COSZ_TEST .and. type == COSZ_FILL )then
                  status = bundle_closeto_expected( bun, expected_bun, eps=8.e-1_r8, crit_type="rms_diff" )
                  call test_is( status, "Test if resultant bundle is close to what is expected :: "//trim(name)//" "//trim(alogo) )
                  if ( .not. status ) call shr_sys_abort( "bad test" )
               else if ( n == COSZ_TEST .and. type == INDX_FILL )then
                  status = bundle_closeto_expected( bun, expected_bun, eps=1.5e-1_r8, crit_type="rel_diff" )
                  call test_is( status, "Test if resultant bundle is close to what is expected :: "//trim(name)//" "//trim(alogo) )
                  if ( .not. status ) call shr_sys_abort( "bad test" )
               else if ( n == COSZ_TEST )then
                  status = bundle_closeto_expected( bun, expected_bun, eps=5.e+1_r8, crit_type="rms_diff" )
                  call test_is( status, "Test if resultant bundle is close to what is expected :: "//trim(name)//" "//trim(alogo) )
               else
                  status = bundle_is_expected(      bun, expected_bun )
                  call test_is( status, "Test if resultant bundle is what is expected :: "//trim(name)//" "//trim(alogo) )
               end if
               call shr_date_adv1step( sdate )
            end do
         end do
         ! Clean the bundles
         call dshr_bundle_clean( bun        )
         call dshr_bundle_clean( expected_bun )
         call dshr_bundle_clean( bunlb        )
         call dshr_bundle_clean( bunub        )
      end do
   end do
   call dshr_domain_clean( domain )

   !
   ! Now do some tests that should fail
   !
   ! Create different sized domains
   dom_name = "test_domain"
   ni = 100
   nj = 300
   call dshr_domain_fill(domain,dom_name,ni,nj)
   dom_name = "test_domain2"
   ni = 10
   nj = 30
   call dshr_domain_fill(domain2,dom_name,ni,nj)

   call dshr_tInterp_setAbort( .false. )

   filltype = "index"
   do type = 1, num_fail
      if (    type == 1 )then
         alogo    = 'linear'
         name = "bunlb_domain_different"
         fieldList  = "T:U:V:TS:PRECT:FSDS:FLDS"
         fieldList2 = fieldList
         domainP1 => domain2
         domainP2 => domain
         domainP3 => domain
      else if ( type == 2 )then
         alogo    = 'linear'
         name = "bunub_domain_different"
         fieldList  = "T:U:V:TS:PRECT:FSDS:FLDS"
         fieldList2 = fieldList
         domainP1 => domain
         domainP2 => domain2
         domainP3 => domain
      else if ( type == 3 )then
         alogo    = 'linear'
         name = "bun_domain_different"
         fieldList  = "T:U:V:TS:PRECT:FSDS:FLDS"
         fieldList2 = fieldList
         domainP1 => domain
         domainP2 => domain
         domainP3 => domain2
      else if ( type == 4 )then
         alogo    = 'linear'
         name = "ublb_fieldLists_different"
         fieldList  = "T:U:V:TS:PRECT:FSDS:FLDS"
         fieldList2 = "T:U:V:TS:PRECT:FSDS:FLDS:THING"
         domainP1 => domain
         domainP2 => domain
         domainP3 => domain
      else if ( type == 5 )then
         alogo      = 'zztop'
         name = "bad_alogo"
         fieldList  = "T:U:V:TS:PRECT:FSDS:FLDS"
         fieldList2 = fieldList
         domainP1 => domain
         domainP2 => domain
         domainP3 => domain
      else if ( type == 6 )then
         call dshr_bundle_setAbort( .false. )
         alogo    = 'linear'
         name = "dates_not_set"
         fieldList  = "T:U:V:TS:PRECT:FSDS:FLDS"
         fieldList2 = fieldList
         domainP1 => domain
         domainP2 => domain
         domainP3 => domain
      else if ( type == 7 )then
         alogo = 'linear'
         name = "ub date before bundle date"
         fieldList  = "T:U:V:TS:PRECT:FSDS:FLDS"//dshr_tInterp_cosZScaleList
         fieldList2 = fieldList
         domainP1 => domain
         domainP2 => domain
         domainP3 => domain
      else
         call shr_sys_abort( "type is out of bounds for fail tests" )
      end if
      call dshr_bundle_create( bunlb, domainP1,  "lb_"//trim(name), fieldList  )
      call dshr_bundle_create( bunub, domainP2,  "ub_"//trim(name), fieldList2 )
      call dshr_bundle_create( bun,   domainP3,  "rs_"//trim(name), fieldList  )
      date_lb(1) = 19480101
      date_ub(1) = 19490102
      sec_lb(1)  = 0
      sec_ub(1)  = 0
      if (      type /= 6 )then
         call dshr_bundle_putDate( bunlb, cdate=date_lb(1), sec=sec_lb(1) )
         call dshr_bundle_putDate( bunub, cdate=date_ub(1), sec=sec_ub(1) )
         date_in = date_lb(1)
         sec_in  = sec_lb(1)
         call dshr_bundle_putDate( bun,   cdate=date_in, sec=sec_in)
      end if
      if ( type == 7 )then
         call dshr_bundle_putDate( bunlb, cdate=date_ub(1), sec=sec_ub(1) )
         call dshr_bundle_putDate( bunub, cdate=date_lb(1), sec=sec_lb(1) )
         sdateUB = shr_date_initCDate( date_lb(1), num_periods*num_times, sec_lb(1) )
         call dshr_tInterp_calcCosZenNorm( bunlb, sdateUB, orb_eccen,  &
                                           orb_mvelpp, orb_lambm0, orb_obliqr, rc=rc )
         if ( rc == 0 )then
            status = .false.
         else
            status = .true.
         end if
         call test_is( status, "Test that calcCosZenNorm fails when ub date < bundle date: test" )
      end if
      call dshr_bundle_fill( bunlb, type=filltype )
      call dshr_bundle_fill( bunub, type=filltype )
      call dshr_bundle_fill( bun,   type=filltype )
      ! Interpolation
      call dshr_tInterp_data(bunlb,bunub,bun,algo=alogo, rc=rc )
      if ( rc == 0 )then
         status = .false.
      else
         status = .true.
      end if
      call test_is( status, "Test things that should fail actually do: test type="//trim(name) )

      call dshr_bundle_setAbort( .true. )
      ! Clean the bundles and nullify domain pointers
      nullify( domainP1 )
      nullify( domainP2 )
      nullify( domainP3 )
      call dshr_bundle_clean( bun        )
      call dshr_bundle_clean( bunlb        )
      call dshr_bundle_clean( bunub        )
   end do
   call dshr_domain_clean( domain )
   call dshr_domain_clean( domain2 )
   deallocate( lat )
   deallocate( lon )

   call test_final( )

contains

subroutine bundle_linear_alg( bun1, bun2, factors, bunout )

   implicit none

   type(dshr_bundle_bundleType), intent(IN)    :: bun1, bun2
   real(R8)                    , intent(IN)    :: factors(2)
   type(dshr_bundle_bundleType), intent(INOUT) :: bunout

   real(r8), pointer :: dataout(:,:,:), data1(:,:,:), data2(:,:,:)

   call dshr_bundle_assignPtr( bunout, dataout )
   call dshr_bundle_assignPtr( bun1,   data1   )
   call dshr_bundle_assignPtr( bun2,   data2   )

   dataout = data1*factors(1) + data2*factors(2)

   nullify( dataout )
   nullify( data1   )
   nullify( data2   )

end subroutine bundle_linear_alg

subroutine scalebycosz( ni, nj, nf, num_times, dtime, dataIN, lat, lon, calday, &
                        declin, date_start, sec_start, dataOut, f, option )
  implicit none

  integer,  intent(IN) :: ni, nj, nf, num_times, dtime
  real(r8), intent(IN) :: dataIN(ni,nj)
  real(r8), intent(IN) :: lat(ni,nj)
  real(r8), intent(IN) :: lon(ni,nj)
  real(r8), intent(IN) :: calday
  real(r8), intent(IN) :: declin
  integer,  intent(IN) :: date_start, sec_start
  real(r8), intent(INOUT) :: dataOUT(ni,nj,nf)
  integer,  intent(IN) :: f
  character(len=*), intent(IN), optional :: option

  integer :: i, j, k
  real(r8) :: cosz(ni,nj), norm(ni,nj), blcalday, eccf, bldeclin, cosz_prev(ni,nj)
  real(r8) :: res, eps, thresh
  real(r8) :: blcalday2, bldeclin2, zen1, zen2
  integer :: bldate_in, blsec_in
  type(shr_date) :: blsdate
  character(len=256) :: loption

  eps = 0.001_r8
  thresh = 0.01_r8
  loption = "box"
  if ( present(option) ) loption = option
  if ( trim(loption) /= "box" .and. trim(loption) /= "triangle" .and. trim(loption) /= "simpsons" )then
      call shr_sys_abort( "bad option to scalebycosz" )
  end if
  ! Calculate normalization factor
  norm(:,:) = 0.0_r8
  bldate_in = date_start
  blsec_in  = sec_start
  if ( trim(loption) == "triangle" )then
     blsdate = shr_date_initCDate( bldate_in, 3600*24/dtime, blsec_in )
     blcalday = shr_date_getJulian( blsdate, shift=-dtime )
     call shr_orb_decl(blcalday, orb_eccen, orb_mvelpp, orb_lambm0, orb_obliqr, &
                       bldeclin, eccf)
     do j = 1, nj
     do i = 1, ni
        cosz_prev(i,j) = shr_orb_cosz(blcalday,lat(i,j),lon(i,j),bldeclin)
        if ( cosz_prev(i,j) < thresh ) cosz_prev(i,j) = thresh
        if ( cosz_prev(i,j) < eps    ) cosz_prev(i,j) = eps
     end do
     end do
  end if
  if ( abs((num_times*dtime/3600.0_r8) - 6.0_r8) > 0.05 )then
      call shr_sys_abort( "not even close to a 6-hour period" )
  end if
  if ( trim(loption) == "simpsons" )then
     blsdate = shr_date_initCDate( bldate_in, 3600*24*2/dtime, blsec_in )
     do k = 1, num_times
        blcalday = shr_date_getJulian( blsdate )
        call shr_orb_decl(blcalday, orb_eccen, orb_mvelpp, orb_lambm0, orb_obliqr, &
                          bldeclin, eccf)
        do j = 1, nj
        do i = 1, ni
           cosz(i,j) = shr_orb_cosz(blcalday,lat(i,j),lon(i,j),bldeclin)
           if ( cosz(i,j) < thresh ) cosz(i,j) = thresh
           if ( cosz(i,j) < eps    ) cosz(i,j) = eps
           norm(i,j) = norm(i,j) + cosz(i,j)
        end do
        end do
        call shr_date_adv1step( blsdate )
        blcalday = shr_date_getJulian( blsdate )
        call shr_orb_decl(blcalday, orb_eccen, orb_mvelpp, orb_lambm0, orb_obliqr, &
                          bldeclin, eccf)
        do j = 1, nj
        do i = 1, ni
           cosz(i,j) = shr_orb_cosz(blcalday,lat(i,j),lon(i,j),bldeclin)
           if ( cosz(i,j) < thresh ) cosz(i,j) = thresh
           if ( cosz(i,j) < eps    ) cosz(i,j) = eps
           norm(i,j) = norm(i,j) + 4.0_r8*cosz(i,j)
        end do
        end do
        call shr_date_adv1step( blsdate )
        blcalday = shr_date_getJulian( blsdate )
        call shr_orb_decl(blcalday, orb_eccen, orb_mvelpp, orb_lambm0, orb_obliqr, &
                          bldeclin, eccf)
        do j = 1, nj
        do i = 1, ni
           cosz(i,j) = shr_orb_cosz(blcalday,lat(i,j),lon(i,j),bldeclin)
           if ( cosz(i,j) < thresh ) cosz(i,j) = thresh
           if ( cosz(i,j) < eps    ) cosz(i,j) = eps
           norm(i,j) = norm(i,j) + cosz(i,j)
        end do
        end do
     end do

     do j = 1, nj
     do i = 1, ni
        norm(i,j) = norm(i,j) / (6._r8 * real(num_times,r8))
     end do
     end do
  else if ( trim(loption) /= "exact" )then
     blsdate = shr_date_initCDate( bldate_in, 3600*24/dtime, blsec_in )
     do k = 1, num_times + 1
        blcalday = shr_date_getJulian( blsdate )
        call shr_orb_decl(blcalday, orb_eccen, orb_mvelpp, orb_lambm0, orb_obliqr, &
                          bldeclin, eccf)
        do j = 1, nj
        do i = 1, ni
           cosz(i,j) = shr_orb_cosz(blcalday,lat(i,j),lon(i,j),bldeclin)
           if ( cosz(i,j) < thresh ) cosz(i,j) = thresh
           if ( cosz(i,j) < eps    ) cosz(i,j) = 0.0_r8
           if ( trim(loption) == "triangle" )then
             res  = min(cosz(i,j),cosz_prev(i,j)) + 0.5_r8*abs(cosz(i,j)-cosz_prev(i,j))
             cosz(i,j) = res
           end if
           if ( cosz(i,j) > eps )then
              norm(i,j) = norm(i,j) + cosz(i,j)
           end if
           cosz_prev(i,j) = cosz(i,j)
        end do
        end do
        call shr_date_adv1step( blsdate )
     end do

     do j = 1, nj
     do i = 1, ni
        norm(i,j) = norm(i,j) / real(num_times+1,r8)
     end do
     end do
  else
     ! This is bogus
     blsdate = shr_date_initCDate( bldate_in, 24/6, blsec_in )
     blcalday = shr_date_getJulian( blsdate )
     call shr_date_adv1step( blsdate )
     blcalday2 = shr_date_getJulian( blsdate )
     call shr_orb_decl(blcalday2, orb_eccen, orb_mvelpp, orb_lambm0, orb_obliqr, &
                       bldeclin2, eccf)
     do j = 1, nj
     do i = 1, ni
        cosz_prev(i,j) = shr_orb_cosz(blcalday,lat(i,j),lon(i,j),bldeclin)
        zen1 = acos( cosz_prev(i,j) )
        cosz(i,j) = shr_orb_cosz(blcalday2,lat(i,j),lon(i,j),bldeclin2)
        zen2 = acos( cosz(i,j) )
        norm(i,j) = max(0.0_r8,-sin(zen2)) + max(0.0_r8,sin(zen1))
     end do
     end do
  end if

  ! Apply normalization factor and current cosine solar-zenith angle to data
  do j = 1, nj
  do i = 1, ni
     cosz(i,j) = shr_orb_cosz(calday,lat(i,j),lon(i,j),declin)
     if ( cosz(i,j) <= thresh ) cosz(i,j) = thresh
     if ( cosz(i,j) <= eps    ) cosz(i,j) = 0.0_r8
  end do
  end do
  do j = 1, nj
  do i = 1, ni
     if ( cosz(i,j) > eps )then
        dataOUT(i,j,f) = dataIN(i,j)*cosz(i,j)/norm(i,j)
     else
        dataOUT(i,j,f) = 0.0_r8
     end if
  end do
  end do

end subroutine scalebycosz

end program test_dshr_tInterp
