      module mo_sim_dat
      private
      public :: set_sim_dat
      contains
      subroutine set_sim_dat
      use chem_mods, only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass
      use chem_mods, only : diag_map
      use chem_mods, only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use chem_mods, only : pht_alias_lst, pht_alias_mult
      use chem_mods, only : het_lst, extfrc_lst, inv_lst, slvd_lst
      use abortutils, only : endrun
      use mo_tracname, only : solsym
      use chem_mods, only : frc_from_dataset
      use shr_kind_mod,only : r8 => shr_kind_r8
      use cam_logfile, only : iulog
      implicit none
!--------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------
      integer :: ios
      clscnt(:) = (/ 17 , 0 , 0 , 42 , 0 /)
      cls_rxt_cnt(:,1) = (/ 10 , 39 , 0 , 17 /)
      cls_rxt_cnt(:,4) = (/ 22 , 96 , 114 , 42 /)
      solsym(: 59) = (/ 'O3      ','O       ','O1D     ','O2      ','O2_1S   ', &
                        'O2_1D   ','N2O     ','N       ','NO      ','NO2     ', &
                        'NO3     ','HNO3    ','HO2NO2  ','N2O5    ','CH4     ', &
                        'CH3O2   ','CH3OOH  ','CH2O    ','CO      ','H2      ', &
                        'H       ','OH      ','HO2     ','H2O2    ','CLY     ', &
                        'BRY     ','CL      ','CL2     ','CLO     ','OCLO    ', &
                        'CL2O2   ','HCL     ','HOCL    ','CLONO2  ','BRCL    ', &
                        'BR      ','BRO     ','HBR     ','HOBR    ','BRONO2  ', &
                        'CH3CL   ','CH3BR   ','CFC11   ','CFC12   ','CFC113  ', &
                        'HCFC22  ','CCL4    ','CH3CCL3 ','CF3BR   ','CF2CLBR ', &
                        'CO2     ','N2p     ','O2p     ','Np      ','Op      ', &
                        'NOp     ','e       ','N2D     ','H2O     ' /)

      adv_mass(: 59) = (/ 47.99819946_r8, 15.99940014_r8, 15.99940014_r8, 31.99880028_r8, 31.99880028_r8, &
                          31.99880028_r8, 44.01287842_r8, 14.00673962_r8, 30.00613976_r8, 46.00553894_r8, &
                          62.00494003_r8, 63.01234055_r8, 79.01174164_r8, 108.0104828_r8, 16.04059982_r8, &
                          47.03200150_r8, 48.03939819_r8, 30.02519989_r8, 28.01040077_r8, 2.014800072_r8, &
                          1.007400036_r8, 17.00679970_r8, 33.00619888_r8, 34.01359940_r8, 100.9168472_r8, &
                          99.71685028_r8, 35.45270157_r8, 70.90540314_r8, 51.45209885_r8, 67.45149994_r8, &
                          102.9041977_r8, 36.46009827_r8, 52.45949936_r8, 97.45764160_r8, 115.3566971_r8, &
                          79.90399933_r8, 95.90339661_r8, 80.91139984_r8, 96.91079712_r8, 141.9089355_r8, &
                          50.48590088_r8, 94.93720245_r8, 137.3675079_r8, 120.9132080_r8, 187.3753052_r8, &
                          86.46790314_r8, 153.8217926_r8, 133.4022980_r8, 148.9102020_r8, 165.3645020_r8, &
                          44.00979996_r8, 28.01347923_r8, 31.99880028_r8, 14.00673962_r8, 15.99940014_r8, &
                          30.00613976_r8, 0.5485669826E-03_r8, 14.00673962_r8, 18.01420021_r8 /)

      fix_mass(: 2) = (/ 0.00000000_r8, 28.0134792_r8 /)
      clsmap(: 17,1) = (/ 15, 7, 19, 20, 41, 42, 43, 44, 45, 46, &
                            47, 48, 49, 50, 51, 25, 26 /)
      clsmap(: 42,4) = (/ 1, 2, 3, 4, 5, 6, 8, 9, 10, 22, &
                            11, 12, 13, 14, 16, 17, 18, 21, 23, 24, &
                            59, 27, 28, 29, 30, 31, 32, 33, 34, 35, &
                            36, 37, 38, 39, 40, 52, 53, 54, 55, 56, &
                            58, 57 /)
      permute(: 42,4) = (/ 42, 32, 36, 29, 3, 2, 22, 33, 35, 38, &
                             39, 17, 9, 7, 23, 8, 27, 28, 37, 14, &
                             26, 40, 5, 31, 4, 1, 41, 24, 25, 6, &
                             34, 30, 16, 18, 13, 15, 19, 10, 11, 20, &
                             12, 21 /)
      diag_map(: 42) = (/ 1, 4, 7, 9, 12, 14, 17, 23, 29, 36, &
                            43, 49, 54, 62, 71, 80, 86, 92, 100, 109, &
                           120, 131, 138, 148, 158, 167, 175, 187, 205, 221, &
                           244, 279, 302, 321, 345, 364, 389, 415, 433, 455, &
                           478, 496 /)
      het_lst(: 59) = (/ 'O3      ','O       ','O1D     ','O2      ','O2_1S   ', &
                         'O2_1D   ','N2O     ','N       ','NO      ','NO2     ', &
                         'NO3     ','HNO3    ','HO2NO2  ','N2O5    ','CH4     ', &
                         'CH3O2   ','CH3OOH  ','CH2O    ','CO      ','H2      ', &
                         'H       ','OH      ','HO2     ','H2O2    ','CLY     ', &
                         'BRY     ','CL      ','CL2     ','CLO     ','OCLO    ', &
                         'CL2O2   ','HCL     ','HOCL    ','CLONO2  ','BRCL    ', &
                         'BR      ','BRO     ','HBR     ','HOBR    ','BRONO2  ', &
                         'CH3CL   ','CH3BR   ','CFC11   ','CFC12   ','CFC113  ', &
                         'HCFC22  ','CCL4    ','CH3CCL3 ','CF3BR   ','CF2CLBR ', &
                         'CO2     ','N2p     ','O2p     ','Np      ','Op      ', &
                         'NOp     ','e       ','N2D     ','H2O     ' /)
      extfrc_lst(: 11) = (/ 'NO      ','NO2     ','CO      ','Op      ','O2p     ', &
                            'Np      ','N2p     ','N2D     ','N       ','e       ', &
                            'OH      ' /)
      frc_from_dataset(: 11) = (/ .true., .true., .true., .false., .false., &
                                  .false., .false., .false., .false., .false., &
                                  .false. /)
      inv_lst(: 2) = (/ 'M       ', 'N2      ' /)
      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_map; error = ',ios
         call endrun
      end if
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jo2_a           ', 'jo2_b           ', 'jo3_a           ', 'jo3_b           ', &
                                     'jn2o            ', 'jno             ', 'jno_i           ', 'jno2            ', &
                                     'jn2o5_a         ', 'jn2o5_b         ', 'jhno3           ', 'jno3_a          ', &
                                     'jno3_b          ', 'jho2no2_a       ', 'jho2no2_b       ', 'jch3ooh         ', &
                                     'jch2o_a         ', 'jch2o_b         ', 'jh2o_a          ', 'jh2o_b          ', &
                                     'jh2o_c          ', 'jh2o2           ', 'jcl2            ', 'jclo            ', &
                                     'joclo           ', 'jcl2o2          ', 'jhocl           ', 'jhcl            ', &
                                     'jclono2_a       ', 'jclono2_b       ', 'jbrcl           ', 'jbro            ', &
                                     'jhobr           ', 'jbrono2_a       ', 'jbrono2_b       ', 'jch3cl          ', &
                                     'jccl4           ', 'jch3ccl3        ', 'jcfcl3          ', 'jcf2cl2         ', &
                                     'jcfc113         ', 'jhcfc22         ', 'jch3br          ', 'jcf3br          ', &
                                     'jcf2clbr        ', 'jco2            ', 'jch4_a          ', 'jch4_b          ', &
                                     'jeuv_1          ', 'jeuv_2          ', 'jeuv_3          ', 'jeuv_4          ', &
                                     'jeuv_5          ', 'jeuv_6          ', 'jeuv_7          ', 'jeuv_8          ', &
                                     'jeuv_9          ', 'jeuv_10         ', 'jeuv_11         ', 'jeuv_12         ', &
                                     'jeuv_13         ', 'jeuv_14         ', 'jeuv_15         ', 'jeuv_16         ', &
                                     'jeuv_17         ', 'jeuv_18         ', 'jeuv_19         ', 'jeuv_20         ', &
                                     'jeuv_21         ', 'jeuv_22         ', 'jeuv_23         ', 'jeuv_24         ', &
                                     'jeuv_25         ', 'usr_O_O2        ', 'cph1            ', 'usr_O_O         ', &
                                     'cph18           ', 'cph19           ', 'cph20           ', 'cph21           ', &
                                     'ag2             ', 'cph22           ', 'cph23           ', 'cph24           ', &
                                     'ag1             ', 'cph17           ', 'cph16           ', 'cph29           ', &
                                     'cph25           ', 'cph26           ', 'cph27           ', 'cph28           ', &
                                     'cph8            ', 'cph12           ', 'cph13           ', 'tag_NO2_NO3     ', &
                                     'usr_N2O5_M      ', 'usr_HNO3_OH     ', 'tag_NO2_HO2     ', 'usr_HO2NO2_M    ', &
                                     'usr_CO_OH_b     ', 'cph5            ', 'cph7            ', 'cph15           ', &
                                     'cph3            ', 'cph11           ', 'cph14           ', 'cph4            ', &
                                     'cph9            ', 'usr_HO2_HO2     ', 'tag_CLO_CLO     ', 'usr_CL2O2_M     ', &
                                     'het1            ', 'het2            ', 'het3            ', 'het4            ', &
                                     'het5            ', 'het6            ', 'het7            ', 'het8            ', &
                                     'het9            ', 'het10           ', 'het11           ', 'het12           ', &
                                     'het13           ', 'het14           ', 'het15           ', 'het16           ', &
                                     'het17           ', 'ion1            ', 'ion2            ', 'ion3            ', &
                                     'ion4            ', 'ion5            ', 'ion6            ', 'ion7            ', &
                                     'ion8            ', 'ion9            ', 'ion11           ', 'elec1           ', &
                                     'elec2           ', 'elec3           ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, &
                                       11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, &
                                       31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
                                       41, 42, 43, 44, 45, 46, 47, 48, 49, 50, &
                                       51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
                                       61, 62, 63, 64, 65, 66, 67, 68, 69, 70, &
                                       71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
                                       82, 83, 84, 85, 86, 87, 88, 89, 108, 109, &
                                      110, 111, 114, 115, 116, 119, 120, 122, 127, 129, &
                                      138, 139, 140, 142, 144, 145, 146, 151, 152, 153, &
                                      171, 172, 202, 203, 204, 205, 206, 207, 208, 209, &
                                      210, 211, 212, 213, 214, 215, 216, 217, 218, 219, &
                                      220, 221, 222, 223, 224, 225, 226, 227, 229, 230, &
                                      231, 232 /)
      if( allocated( pht_alias_lst ) ) then
         deallocate( pht_alias_lst )
      end if
      allocate( pht_alias_lst(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_lst; error = ',ios
         call endrun
      end if
      if( allocated( pht_alias_mult ) ) then
         deallocate( pht_alias_mult )
      end if
      allocate( pht_alias_mult(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_mult; error = ',ios
         call endrun
      end if
      pht_alias_lst(:,1) = (/ 'userdefined     ', 'userdefined     ', '                ', '                ', &
                              '                ', 'userdefined     ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ' /)
      pht_alias_lst(:,2) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ', 'userdefined     ', 'userdefined     ', 'userdefined     ', &
                              'userdefined     ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8 /)
      end subroutine set_sim_dat
      end module mo_sim_dat
