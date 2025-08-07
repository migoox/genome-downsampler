include_guard(GLOBAL)
include(FetchContent)

macro(source_fetcher_begin name)
  set(EXPORT_COMPILE_COMMANDS_TEMP ${CMAKE_EXPORT_COMPILE_COMMANDS})
  set(BUILD_SHARED_LIBS_TEMP ${BUILD_SHARED_LIBS})

  set(CMAKE_EXPORT_COMPILE_COMMANDS OFF CACHE BOOL "" FORCE)
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

  message(CHECK_START "Fetching ${name} source")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")
endmacro()

macro(source_fetcher_end)
  list(POP_BACK CMAKE_MESSAGE_INDENT)
  message(CHECK_PASS "fetched")

  set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_TEMP} CACHE BOOL "" FORCE)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ${EXPORT_COMPILE_COMMANDS_TEMP} CACHE BOOL "" FORCE)
endmacro()


# =============================================================================
# CLI11 
# =============================================================================
function(fetch_source_CLI11)
    if (NOT TARGET CLI11)
        source_fetcher_begin("CLI11")

        FetchContent_Declare(
            CLI11
            GIT_REPOSITORY "https://github.com/CLIUtils/CLI11.git"
            GIT_TAG "v2.5.0"
            #   PATCH_COMMAND git apply --ignore-whitespace
            #                 "${CMAKE_CURRENT_LIST_DIR}/../../patches/CLI11.patch"
            UPDATE_DISCONNECTED 1
        )
        FetchContent_MakeAvailable(CLI11)

        source_fetcher_end()
    else()
        message(STATUS "CLI11 is already available, the fetching process has been abandoned")
    endif()
endfunction()

# =============================================================================
# GoogleTest
# =============================================================================
function(fetch_source_googletest)
  if(NOT TARGET gtest)
    source_fetcher_begin("GoogleTest")

    FetchContent_Declare(
      googletest
      GIT_REPOSITORY "https://github.com/google/googletest.git"
      GIT_TAG "v1.15.2"
    )

    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    source_fetcher_end()
  else()
    message(STATUS "GoogleTest is already available, the fetching process has been abandoned")
  endif()
endfunction()
