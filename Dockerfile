# syntax=docker/dockerfile:1

FROM alpine:latest

# TODO should specify versions--also maybe should have metis source code w/ meshbuilder, can we do that?
RUN apk update && \
    apk add --no-cache \
        build-base \
        cmake \
        bash \
        zsh\
        perl\
        git

WORKDIR /meshbuild
COPY CMakeLists.txt ./
COPY src/ ./src/
COPY meshbuild_workflow.sh .


WORKDIR /meshbuild/build

RUN cmake -B . -S .. && \
    cmake --build . --target all

WORKDIR /meshbuild/src/metis_builds/GKlib

RUN make config prefix=~/src/metis_builds \
    && make install

WORKDIR  /meshbuild/src/metis_builds/METIS
RUN make config prefix=~/src/metis_builds \
    && make install

WORKDIR /meshbuild

# for volume and where meshbuilder workflow lives
RUN mkdir data && mv src/workflow/* data/ && mv build/MeshBuilder data/ \
    && mv src/metis_builds/METIS/build/programs/gpmetis data/




