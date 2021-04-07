#!/bin/sh

make DSP;
./DSP
python3 plot_maps.py
