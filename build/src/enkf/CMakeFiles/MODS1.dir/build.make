# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build

# Include any dependencies generated for this target.
include src/enkf/CMakeFiles/MODS1.dir/depend.make

# Include the progress variables for this target.
include src/enkf/CMakeFiles/MODS1.dir/progress.make

# Include the compile flags for this target's objects.
include src/enkf/CMakeFiles/MODS1.dir/flags.make

src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o: src/enkf/CMakeFiles/MODS1.dir/flags.make
src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o: ../src/enkf/gridinfo_fv3reg.f90
	$(CMAKE_COMMAND) -E cmake_progress_report /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o"
	cd /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build/src/enkf && /apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/bin/intel64/ifort  $(Fortran_DEFINES) $(Fortran_FLAGS) -O3  -warn all -implicitnone -traceback -fp-model strict -convert big_endian -DGFS -D_REAL8_   -DRR_CLOUDANALYSIS -c /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/src/enkf/gridinfo_fv3reg.f90 -o CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o

src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o.requires:
.PHONY : src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o.requires

src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o.provides: src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o.requires
	$(MAKE) -f src/enkf/CMakeFiles/MODS1.dir/build.make src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o.provides.build
.PHONY : src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o.provides

src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o.provides.build: src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o

MODS1: src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o
MODS1: src/enkf/CMakeFiles/MODS1.dir/build.make
.PHONY : MODS1

# Rule to build all files generated by this target.
src/enkf/CMakeFiles/MODS1.dir/build: MODS1
.PHONY : src/enkf/CMakeFiles/MODS1.dir/build

src/enkf/CMakeFiles/MODS1.dir/requires: src/enkf/CMakeFiles/MODS1.dir/gridinfo_fv3reg.f90.o.requires
.PHONY : src/enkf/CMakeFiles/MODS1.dir/requires

src/enkf/CMakeFiles/MODS1.dir/clean:
	cd /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build/src/enkf && $(CMAKE_COMMAND) -P CMakeFiles/MODS1.dir/cmake_clean.cmake
.PHONY : src/enkf/CMakeFiles/MODS1.dir/clean

src/enkf/CMakeFiles/MODS1.dir/depend:
	cd /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/src/enkf /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build/src/enkf /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build/src/enkf/CMakeFiles/MODS1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/enkf/CMakeFiles/MODS1.dir/depend

