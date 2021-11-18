FROM ubuntu:20.04 as builder

RUN apt update && apt install -y --no-install-recommends \
    cmake \
    g++ \
    make \
    libboost-python-dev \
    libboost-numpy-dev \
    python3-numpy \
    python3-pip && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tmp
COPY . .

# Generate python wheel (this calls cmake build internally)
RUN python3 setup.py bdist_wheel && ls dist/

# Using multi-stage build, shrinking final image from 563MB to 151MB
FROM ubuntu:20.04
RUN apt update && apt install -y --no-install-recommends \
    libboost-python1.71.0 \
    libboost-numpy1.71.0 \
    python3-pip \
    python3-numpy && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

# Get python version specific wheel: (cp38 := CPython 3.8)
COPY --from=builder tmp/dist/*cp38*.whl ./
RUN python3 --version && pip install *.whl

COPY --from=builder /tmp/tests/test.py ./
CMD ["python3", "test.py"]

LABEL Name=demag_nonequi_layers Version=0.0.1
