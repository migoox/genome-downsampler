# From: https://raw.githubusercontent.com/genome/build-common/master/cmake/FindHTSlib.cmake
# - Try to find HTSlib
# Once done, this will define
#
#  HTSlib_FOUND - system has HTSlib
#  HTSlib_INCLUDE_DIRS - the HTSlib include directories
#  HTSlib_LIBRARIES - link these to use HTSlib
#  HTSlib::HTSlib - CMake target

set(HTSLIB_SEARCH_DIRS
    ${HTSLIB_SEARCH_DIRS}
    $ENV{HTLSIB_ROOT}
    /gsc/pkg/bio/htslib
    /usr
    /usr/local)

set(_htslib_ver_path "htslib-${htslib_FIND_VERSION}")
include(LibFindMacros)

# Dependencies
libfind_package(HTSlib ZLIB)

# Include dir
find_path(
    HTSlib_INCLUDE_DIR
    NAMES ${HTSLIB_ADDITIONAL_HEADERS} sam.h
    PATHS ${HTSLIB_SEARCH_DIRS}
    PATH_SUFFIXES include include/htslib htslib/${_htslib_ver_path}/htslib
    HINTS ENV HTSLIB_ROOT)

# Finally the library itself
find_library(
    HTSlib_LIBRARY
    NAMES hts libhts.a hts.a
    PATHS ${HTSlib_INCLUDE_DIR} ${HTSLIB_SEARCH_DIRS}
    NO_DEFAULT_PATH
    PATH_SUFFIXES lib lib64 ${_htslib_ver_path}
    HINTS ENV HTSLIB_ROOT)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this lib depends on.
set(HTSlib_PROCESS_INCLUDES HTSlib_INCLUDE_DIR ZLIB_INCLUDE_DIR)
set(HTSlib_PROCESS_LIBS HTSlib_LIBRARY ZLIB_LIBRARIES)
libfind_process(HTSlib)

# Filter the libraries (they should end with .so)
foreach(lib ${HTSlib_LIBRARIES})
    if(${lib} MATCHES "\\.so$")
        list(APPEND filtered_HTSLib_LIBRARIES ${lib})
    endif()
endforeach()

set(HTSlib_LIBRARIES ${filtered_HTSLib_LIBRARIES})

# If path ends with '/include/htslib', use '/include' instead
if(HTSlib_INCLUDE_DIR MATCHES ".*/include/htslib$")
    get_filename_component(HTSlib_INCLUDE_DIR "${HTSlib_INCLUDE_DIR}" DIRECTORY)
endif()

if(HTSlib_FOUND)
    set(HTSlib_INCLUDE_DIRS ${HTSlib_INCLUDE_DIR})
    set(HTSlib_LIBRARIES ${HTSlib_LIBRARY})

    if(NOT TARGET HTSlib::HTSlib)
        add_library(HTSlib::HTSlib UNKNOWN IMPORTED)
        set_target_properties(HTSlib::HTSlib PROPERTIES
            IMPORTED_LOCATION "${HTSlib_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${HTSlib_INCLUDE_DIR}"
        )
    endif()
endif()
