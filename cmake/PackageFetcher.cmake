macro(fetch_package target_name repo tag)
    if (NOT TARGET ${target_name} AND NOT TARGET ${target_name}::${target_name})
        find_package(${target_name} QUIET CONFIG)
        if (NOT ${target_name}_FOUND AND NOT TARGET ${target_name}::${target_name})
            message(STATUS "[fetch_package] ${target_name} not found → downloading ${tag}")

            FetchContent_Declare(
                    ${target_name}
                    GIT_REPOSITORY ${repo}
                    GIT_TAG ${tag}
                    GIT_SHALLOW TRUE
                    FIND_PACKAGE_ARGS CONFIG
            )

            # special handling for VTK
	if (${target_name} STREQUAL "VTK")
                message(STATUS "[fetch_package] Configuring headless VTK build (StandAlone group)...")

                # Disable extras
                set(VTK_BUILD_TESTING OFF CACHE BOOL "" FORCE)
                set(VTK_WRAP_PYTHON OFF CACHE BOOL "" FORCE)
                set(VTK_WRAP_JAVA OFF CACHE BOOL "" FORCE)
                set(VTK_USE_MPI OFF CACHE BOOL "" FORCE) # Change to ON later if your CFD uses parallel MPI
                set(VTK_USE_TK OFF CACHE BOOL "" FORCE)
                set(VTK_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)

                # Disable UI, Rendering, and heavy groups to keep the build fast
                set(VTK_GROUP_ENABLE_Imaging DONT_WANT CACHE STRING "" FORCE)
                set(VTK_GROUP_ENABLE_MPI DONT_WANT CACHE STRING "" FORCE)
                set(VTK_GROUP_ENABLE_Qt DONT_WANT CACHE STRING "" FORCE)
                set(VTK_GROUP_ENABLE_Rendering DONT_WANT CACHE STRING "" FORCE)
                set(VTK_GROUP_ENABLE_Views DONT_WANT CACHE STRING "" FORCE)
                set(VTK_GROUP_ENABLE_Web DONT_WANT CACHE STRING "" FORCE)

                # ENABLE the StandAlone group (This includes Common, Filters, and IO)
                set(VTK_GROUP_ENABLE_StandAlone YES CACHE STRING "" FORCE)
	endif ()
	    	FetchContent_MakeAvailable(${target_name})
        else ()
            message(STATUS "[fetch_package] Found ${target_name}")
        endif ()
    endif ()
endmacro()

include(FetchContent)

function(fetch_gmsh)
    # 1. Prevent re-running if target already exists
    if(TARGET gmsh::gmsh)
        return()
    endif()

    message(STATUS "[fetch_gmsh] Searching for Gmsh...")

    # 2. Try modern CMake Config mode (works if gmshConfig.cmake is in path)
    find_package(gmsh QUIET CONFIG)
    if(gmsh_FOUND AND TARGET gmsh::gmsh)
        message(STATUS "[fetch_gmsh] Found system Gmsh via Config")
        return()
    endif()

    # 3. Manual search for headers/libs (for standard Linux installs missing config files)
    find_path(GMSH_INC NAMES gmsh.h HINTS /usr/include /usr/local/include PATH_SUFFIXES gmsh)
    find_library(GMSH_LIB NAMES gmsh HINTS /usr/lib /usr/lib64 /usr/lib/x86_64-linux-gnu)

    if(GMSH_INC AND GMSH_LIB)
        message(STATUS "[fetch_gmsh] Found system Gmsh (Manual Search)")
        message(STATUS "  - Headers: ${GMSH_INC}")
        message(STATUS "  - Library: ${GMSH_LIB}")
        
        add_library(gmsh::gmsh SHARED IMPORTED GLOBAL)
        set_target_properties(gmsh::gmsh PROPERTIES
            IMPORTED_LOCATION "${GMSH_LIB}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INC}"
        )
        return()
    endif()

    # 4. Fallback: Download official SDK (Safe for servers without sudo)
    message(STATUS "[fetch_gmsh] System files not found. Downloading Gmsh SDK 4.15.1...")

    FetchContent_Declare(
        gmsh_sdk
        URL "file:///home/vva017/gmsh-4.15.1-Linux64-sdk.tgz"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )

    FetchContent_GetProperties(gmsh_sdk)
    if(NOT gmsh_sdk_POPULATED)
        FetchContent_Populate(gmsh_sdk)
    endif()

    # Define paths within the downloaded SDK
    set(GMSH_SDK_ROOT "${gmsh_sdk_SOURCE_DIR}")
    set(GMSH_SDK_INC  "${GMSH_SDK_ROOT}/include")
    set(GMSH_SDK_LIB  "${GMSH_SDK_ROOT}/lib/libgmsh.so")

    # Final verification of the download
    if(EXISTS "${GMSH_SDK_INC}/gmsh.h" AND EXISTS "${GMSH_SDK_LIB}")
        add_library(gmsh::gmsh SHARED IMPORTED GLOBAL)
        set_target_properties(gmsh::gmsh PROPERTIES
            IMPORTED_LOCATION "${GMSH_SDK_LIB}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMSH_SDK_INC}"
        )
        message(STATUS "[fetch_gmsh] Successfully configured Gmsh from SDK")
    else()
        message(FATAL_ERROR "[fetch_gmsh] Failed to find gmsh.h or libgmsh.so in downloaded SDK!")
    endif()
endfunction()

function(fetch_metis)
    if(TARGET metis::metis)
        message(STATUS "[fetch_metis] METIS target already exists")
        return()
    endif()

    message(STATUS "[fetch_metis] Manually fetching and patching METIS...")

    # 1. Fetch GKlib manually (METIS requires it)
    FetchContent_Declare(
        GKlib
        GIT_REPOSITORY https://github.com/KarypisLab/GKlib.git
        GIT_TAG master
    )
    FetchContent_GetProperties(GKlib)
    if(NOT gklib_POPULATED)
        FetchContent_Populate(GKlib)
        add_subdirectory(${gklib_SOURCE_DIR} ${gklib_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()

    # 2. Fetch METIS manually
    FetchContent_Declare(
        metis
        GIT_REPOSITORY https://github.com/KarypisLab/METIS.git
        GIT_TAG v5.2.1
    )
    FetchContent_GetProperties(metis)
    if(NOT metis_POPULATED)
        FetchContent_Populate(metis)

        # --- THE MANUAL PATCHING ---

        # Tell METIS where we just downloaded GKlib
        set(GKLIB_PATH "${gklib_SOURCE_DIR}" CACHE PATH "Path to GKlib" FORCE)

        # Build shared library (recommended)
        set(SHARED ON CACHE BOOL "Build shared library" FORCE)

        # Create the missing directory that "make config" normally makes
        set(METIS_XINC "${metis_SOURCE_DIR}/build/xinclude")
        file(MAKE_DIRECTORY "${METIS_XINC}")

        # Create the configuration header (defaulting to 32-bit integers)
        file(WRITE "${METIS_XINC}/metis.h" "#define IDXTYPEWIDTH 32\n#define REALTYPEWIDTH 32\n")

        # Append the original metis.h content below our configuration
        file(READ "${metis_SOURCE_DIR}/include/metis.h" METIS_ORIG_H)
        file(APPEND "${METIS_XINC}/metis.h" "${METIS_ORIG_H}")

        # Copy the CMakeLists.txt that METIS explicitly expects to find in xinclude
        file(COPY "${metis_SOURCE_DIR}/include/CMakeLists.txt" DESTINATION "${METIS_XINC}")

        # Finally, add the patched METIS directory to the CMake build
        add_subdirectory(${metis_SOURCE_DIR} ${metis_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()

    if(TARGET metis AND NOT TARGET metis::metis)
   
	target_link_libraries(metis PUBLIC GKlib)       
	add_library(metis::metis ALIAS metis)
        
        # 1. Put build/xinclude FIRST so it finds the configured headers
        target_include_directories(metis PUBLIC
            "${metis_SOURCE_DIR}/build/xinclude"
            "${metis_SOURCE_DIR}/include"
        )
        
        # 2. Force the 32-bit definitions globally to prevent the "idx_t" crash
        target_compile_definitions(metis PUBLIC
            IDXTYPEWIDTH=32
            REALTYPEWIDTH=32
        )
    endif()
    message(STATUS "[fetch_metis] METIS is ready!")
endfunction()
