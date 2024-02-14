cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 17)

set(docker ON) #for docker build

if(NOT docker)
    set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-13) # need to set to gcc compiler to acces _List_iterator
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

set(exe "MeshBuilder") # this line is not required but can be changed to change name of executable
project(${exe})
add_executable(${exe} src/main.cpp)

if(NOT docker)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-ld_classic ")
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fsanitize-address-use-after-scope -g") # different options available for using Sanitizer

# currently only tested on new macs (M1 and M2 chips)
message(WARNING "Compiling MeshBuilder for tRIBS")


target_compile_definitions(${exe}
            PRIVATE
            LINUX_32
    )

target_sources(
        ${exe} PRIVATE src/main.cpp
        src/Headers/globalFns.h
        src/Headers/globalFns.cpp
        src/Headers/Definitions.h
        src/Headers/Inclusions.h
        src/Mathutil/mathutil.h
        src/Mathutil/mathutil.cpp
        src/Mathutil/geometry.h
        src/Mathutil/predicates.h
        src/Mathutil/predicates.cpp
        src/bFlowNet/bFlowNet.h
        src/bFlowNet/bFlowNet.cpp
        src/bGraph/bGraph.h
        src/bGraph/bGraph.cpp
        src/bGraph/bGraphNode.h
        src/bGraph/bGraphNode.cpp
        src/bInOut/bInputFile.h
        src/bInOut/bInputFile.cpp
        src/bMeshBuilder/bMeshBuilder.h
        src/bMeshBuilder/bMeshBuilder.cpp
        src/bMeshElements/bMeshElements.h
        src/bMeshElements/bMeshElements.cpp
        src/bMesh/bTriangulator.h
        src/bMesh/bTriangulator.cpp
#
        src/bMesh/heapsort.h
        src/bMeshList/bMeshList.h
        src/bMesh/bMesh.h

)