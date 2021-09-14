#!/bin/bash

cd $( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # make sure we are in correct dir

mkdir -p build
[ "$1" == "-r" ] && rm -r build/*

(cd build && cmake ../demag_nonequidistant && make -j)

PYTHONPATH=build/:${PYTHONPATH} python3 test.py
