FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

RUN echo "Installing dependencies..." && \
      apt-get -y update && \
      apt-get install -y build-essential \
      git \
      wget \
      vim \
      cmake \
      gdb \
      tmux \
      libboost-all-dev \
      python3-setuptools \
      python3-dev \
      python3-pip

RUN python3 -m pip install tensorflow numpy networkx plotly

WORKDIR /home

RUN git clone -b 'v0.2.14-p0' --single-branch --depth 1 https://github.com/Dobiasd/FunctionalPlus && \
cd FunctionalPlus && \
mkdir -p build && cd build && \
cmake .. && \
make && make install && \
cd ../..

RUN git clone -b '3.3.9' --single-branch --depth 1 https://gitlab.com/libeigen/eigen.git && \
cd eigen && \
mkdir -p build && cd build && \
cmake .. && \
make && make install && \
ln -s /usr/local/include/eigen3/Eigen /usr/local/include/Eigen && \
cd ../..

RUN git clone -b 'v3.9.1' --single-branch --depth 1 https://github.com/nlohmann/json && \
cd json && \
mkdir -p build && cd build && \
cmake -DBUILD_TESTING=OFF .. && \
make && make install && \
cd ../..

RUN git clone https://github.com/Dobiasd/frugally-deep && \
cd frugally-deep && \
mkdir -p build && cd build && \
cmake .. && \
make && make install && \
cd ../..

ENV CUDA_VISIBLE_DEVICES=''

