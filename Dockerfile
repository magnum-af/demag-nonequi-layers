# GCC support can be specified at major, minor, or micro version
# (e.g. 8, 8.2 or 8.2.0).
# See https://hub.docker.com/r/library/gcc/ for all supported GCC
# tags from Docker Hub.
# See https://docs.docker.com/samples/library/gcc/ for more on how to use this image
# FROM gcc:latest
FROM ubuntu:latest

RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y \
    cmake \
    g++ \
    libboost-python-dev \
    libboost-numpy-dev \
    python3-pip && \
    pip3 install numpy

# These commands copy your files into the specified directory in the image
# and set that as the working location
COPY . /usr/src/myapp
WORKDIR /usr/src/myapp

# This command compiles your app using GCC, adjust for your source code
# RUN ./run.sh -r
RUN mkdir build && (cd build && cmake .. && make -j)
# RUN g++ -o myapp main.cpp

# This command runs your application, comment out this line to compile only
CMD PYTHONPATH=build python3 test.py
# CMD ["PYTHONPATH=build python3 test.py"]
#CMD ["./test.py"]
#CMD ["./myapp"]

LABEL Name=demagnonequidistant Version=0.0.1
