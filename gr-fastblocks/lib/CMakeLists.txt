# Copyright 2011,2012 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.

########################################################################
# Setup library
########################################################################
include(GrPlatform) #define LIB_SUFFIX

include_directories(${Boost_INCLUDE_DIR}
					${CMAKE_CURRENT_BINARY_DIR}/../include
					${GR_INCLUDE_DIR}
					/home/abhishek/gnuradio/volk/include
					/home/abhishek/gnuradio/build/volk/include
)

link_directories(${Boost_LIBRARY_DIRS})

list(APPEND fastblocks_sources
    ofdm_metric_impl.cc
)

set(fastblocks_sources "${fastblocks_sources}" PARENT_SCOPE)
if(NOT fastblocks_sources)
	MESSAGE(STATUS "No C++ sources... skipping lib/")
	return()
endif(NOT fastblocks_sources)

add_library(gnuradio-fastblocks SHARED ${fastblocks_sources})
target_link_libraries(gnuradio-fastblocks ${Boost_LIBRARIES} ${GNURADIO_ALL_LIBRARIES} /home/abhishek/gnuradio/build/volk/lib/libvolk.so)
set_target_properties(gnuradio-fastblocks PROPERTIES DEFINE_SYMBOL "gnuradio_fastblocks_EXPORTS")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O3 -funroll-loops -msse2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -funroll-loops -msse2")

if(APPLE)
    set_target_properties(gnuradio-fastblocks PROPERTIES
        INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
    )
endif(APPLE)

########################################################################
# Install built library files
########################################################################
install(TARGETS gnuradio-fastblocks
    LIBRARY DESTINATION lib${LIB_SUFFIX} # .so/.dylib file
    ARCHIVE DESTINATION lib${LIB_SUFFIX} # .lib file
    RUNTIME DESTINATION bin              # .dll file
)

########################################################################
# Build and register unit test
########################################################################
include(GrTest)

include_directories(${CPPUNIT_INCLUDE_DIRS})

list(APPEND test_fastblocks_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/test_fastblocks.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qa_fastblocks.cc
)

add_executable(test-fastblocks ${test_fastblocks_sources})

target_link_libraries(
  test-fastblocks
  ${GNURADIO_RUNTIME_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CPPUNIT_LIBRARIES}
  ${GNURADIO_ALL_LIBRARIES}
  gnuradio-fastblocks
  /home/abhishek/gnuradio/build/volk/lib/libvolk.so
)

GR_ADD_TEST(test_fastblocks test-fastblocks)
