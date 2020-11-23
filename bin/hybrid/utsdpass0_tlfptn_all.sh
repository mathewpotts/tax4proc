#!/bin/bash

INPUT=$HOME/data/sdtax4/utsdpass0/listfiles/*
OUTPUT=$HOME/data/sdtax4/utsdpass0/dst/

# utsdpass0 loop
for f in $INPUT
do
    echo "Processing $f..."
    utsdpass0.run -f -incomplete -i $f -o $OUTPUT -d ${f:54:6} > $OUTPUT/SD${f:54:6}.log
done

rm *.Y

# Generate list of files for tlfptn.run
readlink -f $OUTPUT/*.sdpass0.dst.gz > $HOME/data/sdtax4/utsdpass0/want_tlfptn_all.txt

# tlfptn 
tlfptn.run -f -i want_tlfptn_all.txt -o $OUTPUT -det TAX4SD_208 > $HOME/data/sdtax4/utsdpass0/tlfptn_all.log
    
