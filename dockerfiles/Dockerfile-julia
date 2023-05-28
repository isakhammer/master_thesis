FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV HOME_DIR /root

# Essentials
RUN apt-get update && apt-get install -yq \
    apt-utils \
    dialog \
    apt-transport-https \
    wget \
    libasound2 \
    build-essential \
    software-properties-common \
    curl && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get -y upgrade

## Install Juliaup
ENV JULIA_VERSION=1.8.5

# Installation procedure from https://github.com/JuliaLang/juliaup
RUN curl -fsSL https://install.julialang.org -o install_julia.sh \
    && sh install_julia.sh --yes
ENV PATH=/root/.juliaup/bin:$PATH
RUN . /root/.bashrc \
    && cat /root/.bashrc \
    && juliaup add $JULIA_VERSION  \
    && juliaup default $JULIA_VERSION

# Set the working directory
RUN mkdir -p $HOME_DIR/code
ENV WORK_DIR $HOME_DIR/code
WORKDIR $WORK_DIR

# Copy the current directory contents into the container
COPY . .

ENV JULIA_DIR $WORK_DIR/master_thesis/julia
WORKDIR $JULIA_DIR

# Install necessary packages and run the script
RUN julia --project=@. -e 'using Pkg; Pkg.instantiate(); Pkg.resolve()'
RUN julia --project=@. optimize=1 biharmonic_CutFEM_experiments/eoc_test.jl

# Set default command
# CMD ["julia", "--project=@.", "-e", "include('my_script.jl')"]
