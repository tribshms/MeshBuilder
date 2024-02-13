# syntax=docker/dockerfile:1

FROM ubuntu:latest

# Set the working directory
WORKDIR /meshbuilder_wf

# Copy the necessary files
COPY cmake-build-debug-brew_gcc-13/meshBuilder .
COPY src/metis_scripts/gpmetis .
COPY src/metis_scripts/connectivity2metis.pl .
COPY src/metis_scripts/metis2tribs.pl .
COPY src/metis_scripts/run_metis.zsh .
COPY src/metis_scripts/meshbuild_workflow.sh .

# Set executable permissions
RUN chmod +x meshBuilder gpmetis connectivity2metis.pl metis2tribs.pl run_metis.zsh meshbuild_workflow.sh

# Run the workflow script by default when the container starts
CMD ["./meshbuild_workflow.sh"]



#FROM alpine:3.17.0 AS build
#
#RUN apk update && \
#    apk add --no-cache \
#        build-base=0.5-r3 \
#        cmake=3.24.4-r0
#
## build meshbuilder
## Set the working directory
#WORKDIR /meshbuilder
#
#COPY src/ ./src/src/
#COPY CMakeLists.txt ./src
#
## Create a build directory within the working directory
#RUN mkdir build
#
## Set the working directory to the build directory
#WORKDIR /meshbuilder/build
#
## Run CMake from the source directory
#RUN cmake -B . -S ../src && \
#    cmake --build . --target all
#
## copy METIS & related scripts

