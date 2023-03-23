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

ENV CUDA_VISIBLE_DEVICES=''

