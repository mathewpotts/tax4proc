#!/bin/bash

DET=$1

data_hist.py -det $DET -concat
mc_hist_E3.py -det $DET -concat
aperture.py -det $DET
exposure.py -det $DET
spectrum.py -det $DET