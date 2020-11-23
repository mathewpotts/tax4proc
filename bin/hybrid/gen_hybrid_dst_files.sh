#!/bin/bash

START=$1
END=$2

if [[ $# -eq 0 ]]; then
    echo "Usage: gen_hybrid_dst_files.sh YYMMDD YYMMDD\n"
    exit 1
fi

MAIN=$HOME/data
TAX4SD=$MAIN/sdtax4
HYBRID_MDTAX4=$MAIN/hyb_mdtax4
HYBRID_BRTAX4=$MAIN/hyb_brtax4

UTSDPASS0=$TAX4SD/utsdpass0
UTSDPASS0_LISTFILES=$UTSDPASS0/listfiles/
UTSDPASS0_INPUT=$UTSDPASS0/listfiles/*
UTSDPASS0_OUTPUT=$UTSDPASS0/dst/

FD_BYDAY_DATA=/media/glados/Data_Storage_2/hybrid/fd_byday/

TMATCH=/media/glados/Data_Storage_2/hybrid/tmatch/
TMATCH_DST=$TMATCH/dst/
TMATCH_ROOT=$TMATCH/root/

### prep list files for the SDs ###
#$TAX4SD/prep_listfiles.sh $START $END

### utsdpass0 ###
#$UTSDPASS0/utsdpass0_tlfptn_all.sh

### FD_BYDAY ###
#$FD_BYDAY_DATA/fd_byday.sh

# Perform time matching and tlfptn to all time match files
#$TMATCH/tmatch_all.sh

# Pass5 reconstruction
dev_stpfl12_main -det 37 -db -fit 7 -xm 685,720,755,790,825,860 $TMATCH/tmatch_mo_int2_all.tlhbgeomnp.dst.gz -o $TMATCH/pass5/hybrid_${START}_${END}.tlhbgeomnp.dst.gz > $TMATCH/pass5/hybrid_${START}_${END}.tlhbgeomnp.log
dst2rt_ta.run -f $TMATCH/pass5/hybrid_${START}_${END}.tlhbgeomnp.dst.gz -o $TMATCH/pass5
