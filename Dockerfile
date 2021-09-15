FROM ubuntu:latest

RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y \
    cmake \
    g++ \
    libboost-python-dev \
    libboost-numpy-dev \
    python3-pip && \
    pip3 install numpy

ENV APP_PATH=/usr/src/demag_nonequi
COPY . ${APP_PATH}
WORKDIR ${APP_PATH}

RUN mkdir build && (cd build && cmake .. && make -j)

ENV PYTHONPATH=${APP_PATH}/build:${PYTHONPATH}
CMD ["python3", "test.py"]

LABEL Name=demagnonequidistant Version=0.0.1