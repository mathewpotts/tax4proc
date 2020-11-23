#!/bin/bash

SD_PASS0_DATA=/media/glados/Data_Storage_2/hybrid/utsdpass0/dst/
FD_BYDAY_DATA=/media/glados/Data_Storage_2/hybrid/fd_byday/
TMATCH_DST=/media/glados/Data_Storage_2/hybrid/tmatch/dst/
TMATCH_ROOT=/media/glados/Data_Storage_2/hybrid/tmatch/root/
TMATCH=/media/glados/Data_Storage_2/hybrid/tmatch/

for f in $(readlink -f $FD_BYDAY_DATA/*dst.gz)
do
    
    DATE=${f:47:6}

    # Check to see if an SD file for specific date
    if [ $(ls -l $SD_PASS0_DATA/*$DATE* | wc -l) -eq 0 ]; then
	continue
    fi

    # Output only matched events with a 2 sec window to dst file
    tmatch.run -ist 2 -matched_only -ifl1 $f -ifl2 $SD_PASS0_DATA/SD$DATE.incomplete.sdpass0.dst.gz

    # Analyize waveforms for tmatched files
    tlfptn.run tmatch_$DATE.dst.gz


    dst2rt_ta.run -f tmatch_$DATE.tlfptn.dst.gz -o $TMATCH_ROOT
    
done


rm -rf tmatch_mo_int2_all.tlfptn.dst.gz tmatch_mo_int2_all.tlhbgeomnp.dst.gz
mv -f *dst.gz $TMATCH_DST
dstcat.run $TMATCH_DST/*tlfptn.dst.gz -o $TMATCH/tmatch_mo_int2_all.tlfptn.dst.gz
tlhbgeomnp_main.run -sd TASD_ALL -fd TAX4_MIDDLE_DRUM tmatch_mo_int2_all.tlfptn.dst.gz

dst2rt_ta.run -f tmatch_mo_int2_all.tlfptn.dst.gz -o $TMATCH_ROOT

