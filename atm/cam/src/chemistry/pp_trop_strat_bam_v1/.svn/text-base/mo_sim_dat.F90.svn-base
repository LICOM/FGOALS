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
      clscnt(:) = (/ 19 , 0 , 0 , 114 , 0 /)
      cls_rxt_cnt(:,1) = (/ 48 , 40 , 0 , 19 /)
      cls_rxt_cnt(:,4) = (/ 15 , 108 , 209 , 114 /)
      solsym(:133) = (/ 'O3      ','O       ','O1D     ','N2O     ','N       ', &
                        'NO      ','NO2     ','NO3     ','HNO3    ','HO2NO2  ', &
                        'N2O5    ','CH4     ','CH3O2   ','CH3OOH  ','CH3OH   ', &
                        'CH2O    ','CO      ','H2      ','H       ','OH      ', &
                        'HO2     ','H2O2    ','CLY     ','BRY     ','CL      ', &
                        'CL2     ','CLO     ','OCLO    ','CL2O2   ','HCL     ', &
                        'HOCL    ','CLONO2  ','BRCL    ','BR      ','BRO     ', &
                        'HBR     ','HOBR    ','BRONO2  ','CH3CL   ','CH3BR   ', &
                        'CFC11   ','CFC12   ','CFC113  ','HCFC22  ','CCL4    ', &
                        'CH3CCL3 ','CF3BR   ','CF2CLBR ','H2O     ','C2H5OH  ', &
                        'C2H4    ','EO      ','EO2     ','CH3COOH ','GLYALD  ', &
                        'C2H6    ','C2H5O2  ','C2H5OOH ','CH3CHO  ','CH3CO3  ', &
                        'CH3COOOH','C3H6    ','C3H8    ','C3H7O2  ','C3H7OOH ', &
                        'PO2     ','POOH    ','CH3COCH3','RO2     ','ROOH    ', &
                        'BIGENE  ','ENEO2   ','MEK     ','MEKO2   ','MEKOOH  ', &
                        'BIGALK  ','ALKO2   ','ALKOOH  ','ISOP    ','ISOPO2  ', &
                        'ISOPOOH ','MVK     ','MACR    ','MACRO2  ','MACROOH ', &
                        'MCO3    ','HYDRALD ','HYAC    ','CH3COCHO','XO2     ', &
                        'XOOH    ','C10H16  ','TERPO2  ','TERPOOH ','TOLUENE ', &
                        'CRESOL  ','TOLO2   ','TOLOOH  ','XOH     ','BIGALD  ', &
                        'GLYOXAL ','PAN     ','ONIT    ','MPAN    ','ISOPNO3 ', &
                        'ONITR   ','CB1     ','CB2     ','OC1     ','OC2     ', &
                        'SOA     ','SO2     ','DMS     ','SO4     ','NH3     ', &
                        'NH4     ','NH4NO3  ','SSLT01  ','SSLT02  ','SSLT03  ', &
                        'SSLT04  ','DST01   ','DST02   ','DST03   ','DST04   ', &
                        'Rn      ','Pb      ','CO2     ','HCN     ','CH3CN   ', &
                        'C2H2    ','HCOOH   ','HOCH2OO ' /)
      adv_mass(:133) = (/ 47.99820_r8, 15.99940_r8, 15.99940_r8, 44.01288_r8, 14.00674_r8, &
                          30.00614_r8, 46.00554_r8, 62.00494_r8, 63.01234_r8, 79.01174_r8, &
                          108.0105_r8, 16.04060_r8, 47.03200_r8, 48.03940_r8, 32.04000_r8, &
                          30.02520_r8, 28.01040_r8, 2.014800_r8, 1.007400_r8, 17.00680_r8, &
                          33.00620_r8, 34.01360_r8, 100.9168_r8, 99.71685_r8, 35.45270_r8, &
                          70.90540_r8, 51.45210_r8, 67.45150_r8, 102.9042_r8, 36.46010_r8, &
                          52.45950_r8, 97.45764_r8, 115.3567_r8, 79.90400_r8, 95.90340_r8, &
                          80.91140_r8, 96.91080_r8, 141.9089_r8, 50.48590_r8, 94.93720_r8, &
                          137.3675_r8, 120.9132_r8, 187.3753_r8, 86.46790_r8, 153.8218_r8, &
                          133.4023_r8, 148.9102_r8, 165.3645_r8, 18.01420_r8, 46.06580_r8, &
                          28.05160_r8, 61.05780_r8, 77.05720_r8, 60.05040_r8, 60.05040_r8, &
                          30.06640_r8, 61.05780_r8, 62.06520_r8, 44.05100_r8, 75.04240_r8, &
                          76.04980_r8, 42.07740_r8, 44.09220_r8, 75.08360_r8, 76.09100_r8, &
                          91.08300_r8, 92.09040_r8, 58.07680_r8, 89.06820_r8, 90.07560_r8, &
                          56.10320_r8, 105.1088_r8, 72.10260_r8, 103.0940_r8, 104.1014_r8, &
                          72.14380_r8, 103.1352_r8, 104.1426_r8, 68.11420_r8, 117.1198_r8, &
                          118.1272_r8, 70.08780_r8, 70.08780_r8, 119.0934_r8, 120.1008_r8, &
                          101.0792_r8, 100.1130_r8, 74.07620_r8, 72.06140_r8, 149.1186_r8, &
                          150.1260_r8, 136.2284_r8, 185.2340_r8, 186.2414_r8, 92.13620_r8, &
                          108.1356_r8, 173.1406_r8, 174.1480_r8, 190.1474_r8, 98.09820_r8, &
                          58.03560_r8, 121.0479_r8, 119.0743_r8, 147.0847_r8, 162.1179_r8, &
                          147.1259_r8, 12.01100_r8, 12.01100_r8, 12.01100_r8, 12.01100_r8, &
                          144.1320_r8, 64.06480_r8, 62.13240_r8, 96.06360_r8, 17.02894_r8, &
                          18.03634_r8, 80.04128_r8, 58.44247_r8, 58.44247_r8, 58.44247_r8, &
                          58.44247_r8, 135.0640_r8, 135.0640_r8, 135.0640_r8, 135.0640_r8, &
                          222.0000_r8, 207.2000_r8, 44.00980_r8, 27.02514_r8, 41.05094_r8, &
                          26.03680_r8, 46.02460_r8, 63.03140_r8 /)
      fix_mass(: 3) = (/ 0.00000000_r8, 28.0134792_r8, 31.9988003_r8 /)
      clsmap(: 19,1) = (/ 12, 4, 17, 18, 39, 40, 41, 42, 43, 44, &
                            45, 46, 47, 48, 128, 23, 24, 126, 127 /)
      clsmap(:114,4) = (/ 1, 2, 3, 5, 6, 7, 20, 8, 9, 10, &
                            11, 13, 14, 129, 130, 16, 19, 21, 22, 49, &
                            25, 26, 27, 28, 29, 30, 31, 32, 33, 34, &
                            35, 36, 37, 38, 15, 50, 51, 52, 53, 54, &
                            55, 56, 57, 58, 59, 60, 61, 62, 63, 64, &
                            65, 66, 67, 68, 69, 70, 71, 72, 76, 77, &
                            78, 73, 74, 75, 79, 80, 81, 82, 83, 84, &
                            85, 86, 87, 88, 89, 90, 91, 92, 93, 94, &
                            95, 96, 97, 98, 99, 100, 101, 102, 103, 104, &
                           105, 106, 112, 113, 114, 115, 116, 117, 111, 107, &
                           108, 109, 110, 131, 132, 133, 118, 119, 120, 121, &
                           122, 123, 124, 125 /)
      permute(:114,4) = (/ 105, 108, 76, 32, 106, 110, 112, 103, 81, 51, &
                             37, 101, 52, 33, 21, 114, 84, 107, 70, 111, &
                            113, 29, 104, 27, 17, 102, 83, 85, 34, 93, &
                            109, 61, 78, 64, 66, 36, 43, 30, 47, 68, &
                             67, 22, 74, 38, 87, 100, 57, 86, 23, 80, &
                             45, 77, 62, 75, 88, 46, 18, 39, 19, 73, &
                             71, 53, 72, 40, 79, 98, 63, 96, 92, 97, &
                             41, 99, 44, 91, 94, 95, 42, 69, 90, 54, &
                             25, 26, 60, 48, 31, 58, 59, 55, 49, 65, &
                             89, 82, 24, 35, 1, 20, 2, 3, 4, 5, &
                              6, 7, 8, 28, 50, 56, 9, 10, 11, 12, &
                             13, 14, 15, 16 /)
      diag_map(:114) = (/ 1, 2, 3, 4, 5, 7, 8, 10, 11, 12, &
                            13, 14, 15, 16, 17, 18, 19, 22, 25, 28, &
                            31, 34, 38, 43, 45, 50, 53, 56, 61, 63, &
                            67, 71, 75, 79, 83, 88, 92, 98, 103, 110, &
                           115, 120, 124, 131, 134, 140, 147, 153, 159, 163, &
                           167, 173, 179, 184, 191, 199, 206, 212, 217, 223, &
                           230, 235, 243, 251, 259, 267, 272, 276, 280, 289, &
                           297, 308, 319, 333, 341, 348, 360, 370, 379, 395, &
                           406, 412, 420, 427, 436, 448, 462, 472, 485, 499, &
                           512, 518, 529, 538, 551, 565, 584, 604, 621, 646, &
                           678, 696, 733, 751, 788, 836, 897, 920, 940, 972, &
                           988,1078,1098,1117 /)
      het_lst(:133) = (/ 'O3      ','O       ','O1D     ','N2O     ','N       ', &
                         'NO      ','NO2     ','NO3     ','HNO3    ','HO2NO2  ', &
                         'N2O5    ','CH4     ','CH3O2   ','CH3OOH  ','CH3OH   ', &
                         'CH2O    ','CO      ','H2      ','H       ','OH      ', &
                         'HO2     ','H2O2    ','CLY     ','BRY     ','CL      ', &
                         'CL2     ','CLO     ','OCLO    ','CL2O2   ','HCL     ', &
                         'HOCL    ','CLONO2  ','BRCL    ','BR      ','BRO     ', &
                         'HBR     ','HOBR    ','BRONO2  ','CH3CL   ','CH3BR   ', &
                         'CFC11   ','CFC12   ','CFC113  ','HCFC22  ','CCL4    ', &
                         'CH3CCL3 ','CF3BR   ','CF2CLBR ','H2O     ','C2H5OH  ', &
                         'C2H4    ','EO      ','EO2     ','CH3COOH ','GLYALD  ', &
                         'C2H6    ','C2H5O2  ','C2H5OOH ','CH3CHO  ','CH3CO3  ', &
                         'CH3COOOH','C3H6    ','C3H8    ','C3H7O2  ','C3H7OOH ', &
                         'PO2     ','POOH    ','CH3COCH3','RO2     ','ROOH    ', &
                         'BIGENE  ','ENEO2   ','MEK     ','MEKO2   ','MEKOOH  ', &
                         'BIGALK  ','ALKO2   ','ALKOOH  ','ISOP    ','ISOPO2  ', &
                         'ISOPOOH ','MVK     ','MACR    ','MACRO2  ','MACROOH ', &
                         'MCO3    ','HYDRALD ','HYAC    ','CH3COCHO','XO2     ', &
                         'XOOH    ','C10H16  ','TERPO2  ','TERPOOH ','TOLUENE ', &
                         'CRESOL  ','TOLO2   ','TOLOOH  ','XOH     ','BIGALD  ', &
                         'GLYOXAL ','PAN     ','ONIT    ','MPAN    ','ISOPNO3 ', &
                         'ONITR   ','CB1     ','CB2     ','OC1     ','OC2     ', &
                         'SOA     ','SO2     ','DMS     ','SO4     ','NH3     ', &
                         'NH4     ','NH4NO3  ','SSLT01  ','SSLT02  ','SSLT03  ', &
                         'SSLT04  ','DST01   ','DST02   ','DST03   ','DST04   ', &
                         'Rn      ','Pb      ','CO2     ','HCN     ','CH3CN   ', &
                         'C2H2    ','HCOOH   ','HOCH2OO ' /)
      extfrc_lst(: 5) = (/ 'NO      ','NO2     ','CO      ','SO2     ','CB1     ' /)
      frc_from_dataset(: 5) = (/ .true., .true., .true., .true., .true. /)
      inv_lst(: 3) = (/ 'M       ', 'N2      ', 'O2      ' /)
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
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jo2_b           ', 'jo3_a           ', 'jo3_b           ', 'jn2o            ', &
                                     'jno             ', 'jno2            ', 'jn2o5_a         ', 'jn2o5_b         ', &
                                     'jhno3           ', 'jno3_a          ', 'jno3_b          ', 'jho2no2_a       ', &
                                     'jho2no2_b       ', 'jch3ooh         ', 'jch2o_a         ', 'jch2o_b         ', &
                                     'jch4_a          ', 'jch4_b          ', 'jch3cho         ', 'jpooh           ', &
                                     'jch3co3h        ', 'jpan            ', 'jmpan           ', 'jmacr_a         ', &
                                     'jmacr_b         ', 'jmvk            ', 'jc2h5ooh        ', 'jc3h7ooh        ', &
                                     'jrooh           ', 'jacet           ', 'jmgly           ', 'jxooh           ', &
                                     'jonitr          ', 'jisopooh        ', 'jhyac           ', 'jglyald         ', &
                                     'jmek            ', 'jbigald         ', 'jglyoxal        ', 'jalkooh         ', &
                                     'jmekooh         ', 'jtolooh         ', 'jterpooh        ', 'jh2o_a          ', &
                                     'jh2o_b          ', 'jh2o_c          ', 'jh2o2           ', 'jcl2            ', &
                                     'joclo           ', 'jcl2o2          ', 'jhocl           ', 'jhcl            ', &
                                     'jclono2_a       ', 'jclono2_b       ', 'jbrcl           ', 'jbro            ', &
                                     'jhobr           ', 'jbrono2_a       ', 'jbrono2_b       ', 'jch3cl          ', &
                                     'jccl4           ', 'jch3ccl3        ', 'jcfcl3          ', 'jcf2cl2         ', &
                                     'jcfc113         ', 'jhcfc22         ', 'jch3br          ', 'jcf3br          ', &
                                     'jcf2clbr        ', 'jco2            ', 'usr_O_O2        ', 'usr_O_O         ', &
                                     'o1d_n2          ', 'o1d_o2          ', 'ox_l1           ', 'ox_l2           ', &
                                     'ox_l3           ', 'usr_HO2_HO2     ', 'ox_p1           ', 'tag_NO2_NO3     ', &
                                     'usr_N2O5_M      ', 'tag_NO2_OH      ', 'usr_HNO3_OH     ', 'tag_NO2_HO2     ', &
                                     'usr_HO2NO2_M    ', 'ox_p2           ', 'usr_CO_OH_b     ', 'tag_C2H4_OH     ', &
                                     'ox_l6           ', 'ox_p5           ', 'ox_p4           ', 'tag_CH3CO3_NO2  ', &
                                     'ox_p16          ', 'usr_PAN_M       ', 'tag_C3H6_OH     ', 'ox_l4           ', &
                                     'ox_p9           ', 'ox_p3           ', 'usr_CH3COCH3_OH ', 'ox_p10          ', &
                                     'ox_p15          ', 'ox_l7           ', 'ox_p17          ', 'ox_l8           ', &
                                     'ox_p7           ', 'ox_p8           ', 'usr_MCO3_NO2    ', 'usr_MPAN_M      ', &
                                     'ox_l5           ', 'ox_p6           ', 'soa5            ', 'ox_p14          ', &
                                     'ox_p11          ', 'usr_XOOH_OH     ', 'soa4            ', 'ox_p12          ', &
                                     'soa2            ', 'soa1            ', 'soa3            ', 'ox_p13          ', &
                                     'usr_N2O5_aer    ', 'usr_NO3_aer     ', 'usr_NO2_aer     ', 'usr_SO2_OH      ', &
                                     'usr_DMS_OH      ', 'usr_HO2_aer     ', 'tag_CLO_CLO     ', 'usr_CL2O2_M     ', &
                                     'het1            ', 'het2            ', 'het3            ', 'het4            ', &
                                     'het5            ', 'het6            ', 'het7            ', 'het8            ', &
                                     'het9            ', 'het10           ', 'het11           ', 'het12           ', &
                                     'het13           ', 'het14           ', 'het15           ', 'het16           ', &
                                     'het17           ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, &
                                       11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, &
                                       31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
                                       41, 42, 43, 44, 45, 46, 47, 48, 49, 50, &
                                       51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
                                       61, 62, 63, 64, 65, 66, 67, 68, 69, 70, &
                                       71, 73, 74, 75, 76, 90, 99, 100, 106, 111, &
                                      112, 113, 114, 119, 121, 123, 130, 142, 143, 145, &
                                      152, 153, 158, 164, 166, 167, 169, 174, 177, 178, &
                                      187, 189, 191, 195, 196, 203, 209, 210, 213, 215, &
                                      224, 228, 231, 237, 238, 239, 244, 245, 246, 247, &
                                      251, 252, 253, 255, 257, 261, 286, 287, 317, 318, &
                                      319, 320, 321, 322, 323, 324, 325, 326, 327, 328, &
                                      329, 330, 331, 332, 333 /)
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
      pht_alias_lst(:,1) = (/ 'userdefined     ', '                ', '                ', '                ', &
                              'userdefined     ', '                ', '                ', '                ', &
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
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ' /)
      pht_alias_lst(:,2) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', 'jch3ooh         ', &
                              'jh2o2           ', '                ', 'jpan            ', '                ', &
                              '                ', '                ', 'jch3ooh         ', 'jch3ooh         ', &
                              'jch3ooh         ', '                ', '                ', 'jch3ooh         ', &
                              'jch3cho         ', 'jch3ooh         ', '                ', '                ', &
                              'jacet           ', 'jno2            ', 'jmgly           ', 'jch3ooh         ', &
                              'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ' /)
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
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          .28_r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, .2_r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8 /)
      end subroutine set_sim_dat
      end module mo_sim_dat
