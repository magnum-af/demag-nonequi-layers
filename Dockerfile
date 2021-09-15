FROM ubuntu:latest

RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y \
    cmake \
    g++ \
    libboost-python-dev \
    libboost-numpy-dev \
    python3-pip && \
    pip3 install numpy

WORKDIR /tmp
COPY . .

RUN mkdir build && \
    (cd build && cmake .. && make -j && make install) && \
    rm -r build

ENV PYTHONPATH=/usr/local/lib
CMD ["python3", "test.py"]

LABEL Name=demagnonequidistant Version=0.0.1