# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.22.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.22.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/mailund/Classroom/solutions/gsa/fm-index-c-cmake

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/mailund/Classroom/solutions/gsa/fm-index-c-cmake

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target test
test:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running tests..."
	/usr/local/Cellar/cmake/3.22.2/bin/ctest --force-new-ctest-process $(ARGS)
.PHONY : test

# Special rule for the target test
test/fast: test
.PHONY : test/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/local/Cellar/cmake/3.22.2/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/local/Cellar/cmake/3.22.2/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\" \"library\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components
.PHONY : list_install_components/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/local/Cellar/cmake/3.22.2/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/local/Cellar/cmake/3.22.2/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/local/Cellar/cmake/3.22.2/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/local/Cellar/cmake/3.22.2/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local/fast

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/local/Cellar/cmake/3.22.2/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/local/Cellar/cmake/3.22.2/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/mailund/Classroom/solutions/gsa/fm-index-c-cmake/CMakeFiles /Users/mailund/Classroom/solutions/gsa/fm-index-c-cmake//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/mailund/Classroom/solutions/gsa/fm-index-c-cmake/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named fm

# Build rule for target.
fm: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 fm
.PHONY : fm

# fast build rule for target.
fm/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/build
.PHONY : fm/fast

#=============================================================================
# Target rules for targets named cstr

# Build rule for target.
cstr: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 cstr
.PHONY : cstr

# fast build rule for target.
cstr/fast:
	$(MAKE) $(MAKESILENT) -f src/cstr/CMakeFiles/cstr.dir/build.make src/cstr/CMakeFiles/cstr.dir/build
.PHONY : cstr/fast

#=============================================================================
# Target rules for targets named testlib

# Build rule for target.
testlib: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 testlib
.PHONY : testlib

# fast build rule for target.
testlib/fast:
	$(MAKE) $(MAKESILENT) -f src/test/CMakeFiles/testlib.dir/build.make src/test/CMakeFiles/testlib.dir/build
.PHONY : testlib/fast

#=============================================================================
# Target rules for targets named alphabet_test

# Build rule for target.
alphabet_test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 alphabet_test
.PHONY : alphabet_test

# fast build rule for target.
alphabet_test/fast:
	$(MAKE) $(MAKESILENT) -f src/test/CMakeFiles/alphabet_test.dir/build.make src/test/CMakeFiles/alphabet_test.dir/build
.PHONY : alphabet_test/fast

#=============================================================================
# Target rules for targets named bit_vector_test

# Build rule for target.
bit_vector_test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 bit_vector_test
.PHONY : bit_vector_test

# fast build rule for target.
bit_vector_test/fast:
	$(MAKE) $(MAKESILENT) -f src/test/CMakeFiles/bit_vector_test.dir/build.make src/test/CMakeFiles/bit_vector_test.dir/build
.PHONY : bit_vector_test/fast

#=============================================================================
# Target rules for targets named bwt_test

# Build rule for target.
bwt_test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 bwt_test
.PHONY : bwt_test

# fast build rule for target.
bwt_test/fast:
	$(MAKE) $(MAKESILENT) -f src/test/CMakeFiles/bwt_test.dir/build.make src/test/CMakeFiles/bwt_test.dir/build
.PHONY : bwt_test/fast

#=============================================================================
# Target rules for targets named cstr_test

# Build rule for target.
cstr_test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 cstr_test
.PHONY : cstr_test

# fast build rule for target.
cstr_test/fast:
	$(MAKE) $(MAKESILENT) -f src/test/CMakeFiles/cstr_test.dir/build.make src/test/CMakeFiles/cstr_test.dir/build
.PHONY : cstr_test/fast

#=============================================================================
# Target rules for targets named sa_test

# Build rule for target.
sa_test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 sa_test
.PHONY : sa_test

# fast build rule for target.
sa_test/fast:
	$(MAKE) $(MAKESILENT) -f src/test/CMakeFiles/sa_test.dir/build.make src/test/CMakeFiles/sa_test.dir/build
.PHONY : sa_test/fast

#=============================================================================
# Target rules for targets named sais_test

# Build rule for target.
sais_test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 sais_test
.PHONY : sais_test

# fast build rule for target.
sais_test/fast:
	$(MAKE) $(MAKESILENT) -f src/test/CMakeFiles/sais_test.dir/build.make src/test/CMakeFiles/sais_test.dir/build
.PHONY : sais_test/fast

#=============================================================================
# Target rules for targets named skew_test

# Build rule for target.
skew_test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 skew_test
.PHONY : skew_test

# fast build rule for target.
skew_test/fast:
	$(MAKE) $(MAKESILENT) -f src/test/CMakeFiles/skew_test.dir/build.make src/test/CMakeFiles/skew_test.dir/build
.PHONY : skew_test/fast

src/fasta.o: src/fasta.c.o
.PHONY : src/fasta.o

# target to build an object file
src/fasta.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/fasta.c.o
.PHONY : src/fasta.c.o

src/fasta.i: src/fasta.c.i
.PHONY : src/fasta.i

# target to preprocess a source file
src/fasta.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/fasta.c.i
.PHONY : src/fasta.c.i

src/fasta.s: src/fasta.c.s
.PHONY : src/fasta.s

# target to generate assembly for a file
src/fasta.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/fasta.c.s
.PHONY : src/fasta.c.s

src/fastq.o: src/fastq.c.o
.PHONY : src/fastq.o

# target to build an object file
src/fastq.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/fastq.c.o
.PHONY : src/fastq.c.o

src/fastq.i: src/fastq.c.i
.PHONY : src/fastq.i

# target to preprocess a source file
src/fastq.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/fastq.c.i
.PHONY : src/fastq.c.i

src/fastq.s: src/fastq.c.s
.PHONY : src/fastq.s

# target to generate assembly for a file
src/fastq.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/fastq.c.s
.PHONY : src/fastq.c.s

src/fm.o: src/fm.c.o
.PHONY : src/fm.o

# target to build an object file
src/fm.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/fm.c.o
.PHONY : src/fm.c.o

src/fm.i: src/fm.c.i
.PHONY : src/fm.i

# target to preprocess a source file
src/fm.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/fm.c.i
.PHONY : src/fm.c.i

src/fm.s: src/fm.c.s
.PHONY : src/fm.s

# target to generate assembly for a file
src/fm.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/fm.c.s
.PHONY : src/fm.c.s

src/sam.o: src/sam.c.o
.PHONY : src/sam.o

# target to build an object file
src/sam.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/sam.c.o
.PHONY : src/sam.c.o

src/sam.i: src/sam.c.i
.PHONY : src/sam.i

# target to preprocess a source file
src/sam.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/sam.c.i
.PHONY : src/sam.c.i

src/sam.s: src/sam.c.s
.PHONY : src/sam.s

# target to generate assembly for a file
src/sam.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fm.dir/build.make CMakeFiles/fm.dir/src/sam.c.s
.PHONY : src/sam.c.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... install"
	@echo "... install/local"
	@echo "... install/strip"
	@echo "... list_install_components"
	@echo "... rebuild_cache"
	@echo "... test"
	@echo "... alphabet_test"
	@echo "... bit_vector_test"
	@echo "... bwt_test"
	@echo "... cstr"
	@echo "... cstr_test"
	@echo "... fm"
	@echo "... sa_test"
	@echo "... sais_test"
	@echo "... skew_test"
	@echo "... testlib"
	@echo "... src/fasta.o"
	@echo "... src/fasta.i"
	@echo "... src/fasta.s"
	@echo "... src/fastq.o"
	@echo "... src/fastq.i"
	@echo "... src/fastq.s"
	@echo "... src/fm.o"
	@echo "... src/fm.i"
	@echo "... src/fm.s"
	@echo "... src/sam.o"
	@echo "... src/sam.i"
	@echo "... src/sam.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

