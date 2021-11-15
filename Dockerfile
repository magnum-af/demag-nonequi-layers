FROM ubuntu:20.04 as builder

RUN apt update && apt install -y --no-install-recommends \
    g++ \
    make \
    cmake \
    libboost-python-dev \
    libboost-numpy-dev \
    python3-numpy && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tmp
COPY . .

RUN mkdir build && \
    (cd build && cmake .. && make -j && make install) && \
    rm -r build

# Using multi-stage build, shrinking final image from 555MB to 142MB
# Note: this entire block could be skipped and still works
FROM ubuntu:20.04
RUN apt update && apt install -y --no-install-recommends \
    libboost-python1.71.0 \
    libboost-numpy1.71.0 \
    python3-numpy && \
    rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/local/lib/*.so /usr/local/lib
COPY --from=builder /tmp/*.py ./

ENV PYTHONPATH=/usr/local/lib
CMD ["python3", "test.py"]

LABEL Name=demag_nonequi_layers Version=0.0.1
