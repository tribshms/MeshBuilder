# syntax=docker/dockerfile:1

FROM alpine:3.19.1

# TODO should specify versions--also maybe should have metis source code w/ meshbuilder, can we do that?
RUN apk update && \
    apk add --no-cache \
        build-base=0.5-r3 \
        cmake=3.27.8-r0 \
        bash=5.2.21-r0\
        zsh=5.9-r2\
        perl=5.38.2-r0

WORKDIR /meshbuild
COPY CMakeLists.txt ./
COPY src/ ./src/
COPY meshbuild_workflow.sh .
RUN chmod +x meshbuild_workflow.sh

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



