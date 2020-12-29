#!/usr/bin/env python3

# Import Libs
import ROOT
import TDSTio
import TLPY
import os
import argparse

# PASS5 FILE PATHS FOR EACH EPOCH
epoch0 = os.environ['MDTAX4_DATA_PASS5'] + '/epoch0/'
epoch1 = os.environ['MDTAX4_DATA_PASS5'] + '/epoch1/'
epoch2 = os.environ['MDTAX4_DATA_PASS5'] + '/processing/root_all_files/'
