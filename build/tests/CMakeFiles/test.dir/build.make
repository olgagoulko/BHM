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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/olga/Documents/sampling

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/olga/Documents/sampling/build

# Include any dependencies generated for this target.
include tests/CMakeFiles/test.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/test.dir/flags.make

tests/CMakeFiles/test.dir/basic_unittest.cpp.o: tests/CMakeFiles/test.dir/flags.make
tests/CMakeFiles/test.dir/basic_unittest.cpp.o: ../tests/basic_unittest.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/olga/Documents/sampling/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/CMakeFiles/test.dir/basic_unittest.cpp.o"
	cd /home/olga/Documents/sampling/build/tests && g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test.dir/basic_unittest.cpp.o -c /home/olga/Documents/sampling/tests/basic_unittest.cpp

tests/CMakeFiles/test.dir/basic_unittest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test.dir/basic_unittest.cpp.i"
	cd /home/olga/Documents/sampling/build/tests && g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/olga/Documents/sampling/tests/basic_unittest.cpp > CMakeFiles/test.dir/basic_unittest.cpp.i

tests/CMakeFiles/test.dir/basic_unittest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test.dir/basic_unittest.cpp.s"
	cd /home/olga/Documents/sampling/build/tests && g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/olga/Documents/sampling/tests/basic_unittest.cpp -o CMakeFiles/test.dir/basic_unittest.cpp.s

tests/CMakeFiles/test.dir/basic_unittest.cpp.o.requires:
.PHONY : tests/CMakeFiles/test.dir/basic_unittest.cpp.o.requires

tests/CMakeFiles/test.dir/basic_unittest.cpp.o.provides: tests/CMakeFiles/test.dir/basic_unittest.cpp.o.requires
	$(MAKE) -f tests/CMakeFiles/test.dir/build.make tests/CMakeFiles/test.dir/basic_unittest.cpp.o.provides.build
.PHONY : tests/CMakeFiles/test.dir/basic_unittest.cpp.o.provides

tests/CMakeFiles/test.dir/basic_unittest.cpp.o.provides.build: tests/CMakeFiles/test.dir/basic_unittest.cpp.o

# Object files for target test
test_OBJECTS = \
"CMakeFiles/test.dir/basic_unittest.cpp.o"

# External object files for target test
test_EXTERNAL_OBJECTS =

tests/test: tests/CMakeFiles/test.dir/basic_unittest.cpp.o
tests/test: tests/CMakeFiles/test.dir/build.make
tests/test: src/libsamplinglibrary.a
tests/test: tests/CMakeFiles/test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable test"
	cd /home/olga/Documents/sampling/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/test.dir/build: tests/test
.PHONY : tests/CMakeFiles/test.dir/build

tests/CMakeFiles/test.dir/requires: tests/CMakeFiles/test.dir/basic_unittest.cpp.o.requires
.PHONY : tests/CMakeFiles/test.dir/requires

tests/CMakeFiles/test.dir/clean:
	cd /home/olga/Documents/sampling/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/test.dir/clean

tests/CMakeFiles/test.dir/depend:
	cd /home/olga/Documents/sampling/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/olga/Documents/sampling /home/olga/Documents/sampling/tests /home/olga/Documents/sampling/build /home/olga/Documents/sampling/build/tests /home/olga/Documents/sampling/build/tests/CMakeFiles/test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/test.dir/depend

