# CMake generated Testfile for 
# Source directory: /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF
# Build directory: /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(global_T62 "regression_driver.sh" "global_T62" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(global_T62 PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(global_T62_ozonly "regression_driver.sh" "global_T62_ozonly" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(global_T62_ozonly PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(global_4dvar_T62 "regression_driver.sh" "global_4dvar_T62" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(global_4dvar_T62 PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(global_4denvar_T126 "regression_driver.sh" "global_4denvar_T126" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(global_4denvar_T126 PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(global_fv3_4denvar_T126 "regression_driver.sh" "global_fv3_4denvar_T126" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(global_fv3_4denvar_T126 PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(global_fv3_4denvar_C192 "regression_driver.sh" "global_fv3_4denvar_C192" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(global_fv3_4denvar_C192 PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(global_lanczos_T62 "regression_driver.sh" "global_lanczos_T62" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(global_lanczos_T62 PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(arw_netcdf "regression_driver.sh" "arw_netcdf" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(arw_netcdf PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(
          arw_binary "regression_driver.sh" "
          arw_binary" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(
          arw_binary PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(nmm_binary "regression_driver.sh" "nmm_binary" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(nmm_binary PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(nmm_netcdf "regression_driver.sh" "nmm_netcdf" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(nmm_netcdf PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(nmmb_nems_4denvar "regression_driver.sh" "nmmb_nems_4denvar" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(nmmb_nems_4denvar PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(hwrf_nmm_d2 "regression_driver.sh" "hwrf_nmm_d2" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(hwrf_nmm_d2 PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(hwrf_nmm_d3 "regression_driver.sh" "hwrf_nmm_d3" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(hwrf_nmm_d3 PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(rtma "regression_driver.sh" "rtma" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(rtma PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(global_enkf_T62 "regression_driver.sh" "global_enkf_T62" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(global_enkf_T62 PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(netcdf_fv3_regional "regression_driver.sh" "netcdf_fv3_regional" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(netcdf_fv3_regional PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(global_C96_fv3aero "regression_driver.sh" "global_C96_fv3aero" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(global_C96_fv3aero PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
ADD_TEST(global_C96_fv3aerorad "regression_driver.sh" "global_C96_fv3aerorad" "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build")
SET_TESTS_PROPERTIES(global_C96_fv3aerorad PROPERTIES  TIMEOUT "86400" WORKING_DIRECTORY "/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/regression")
SUBDIRS(src/ncdiag)
SUBDIRS(src/fv3gfs_ncio)
SUBDIRS(src/GSD/gsdcloud)
SUBDIRS(src/gsi)
SUBDIRS(src/enkf)
SUBDIRS(util/ndate)
SUBDIRS(util/EnKF/arw/src)
SUBDIRS(util/Analysis_Utilities/read_diag)
SUBDIRS(util/radar_process/radialwind)
SUBDIRS(util/radar_process/reflectivity)
SUBDIRS(util/bufr_tools)
