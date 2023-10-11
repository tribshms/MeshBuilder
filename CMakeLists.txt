cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-13) # need to set to gcc compiler to acces _List_iterator
set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

set(exe "meshBuilder") # this line is not required but can be changed to change name of executable
project(${exe})
add_executable(${exe} src/main.cpp)


# Set linker flags
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-ld_classic")

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