#!/bin/bash

readlink -f $HOME/data/mdtax4/pass3/data/*/*ps2g.dst > want_tafd_all.ok.txt

ta_byday.run -i want_tafd_all.ok.txt -f -o $HOME/data/hyb_mdtax4/fd_byday/ -v
