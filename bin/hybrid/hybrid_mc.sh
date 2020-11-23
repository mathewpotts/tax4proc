#!/bin/bash

HYBRID_MC=/home/glados/Data_Storage_2/hybrid/hybrid_mc

# Make directories
ls /media/glados/Data_Storage_2/mdtax4/mc/ | grep [0-9][0-9] |awk '{print "mkdir -p "$1"/pass3"}' | bash

# Loop through directories and make SD simulation
for dir in */; do
    cd $dir/pass3

    echo "Adding SD simulation data to $dir FD MC..." 
    readlink -f /media/glados/Data_Storage_2/mdtax4/mc/$dir/pass3/*ps3.dst > want_tamc_all.ok.txt
    tl_1sd_hy_simulation.run -sd_det 1 -i want_tamc_all.ok.txt -f

    echo "Adding tlfptn bank to $dir MC"
    readlink -f /media/glados/Data_Storage_2/hybrid/hybrid_mc/$dir/pass3/*ps3.wsds.dst.gz > want_wsds.txt
    tlfptn.run -det 1 -i want_wsds.txt -f -no_bw

    echo "Running tlhbgeonp on $dir MC"
    readlink -f /media/glados/Data_Storage_2/hybrid/hybrid_mc/$dir/pass3/*ps3.wsds.tlfptn.dst.gz > want_wsds.tlfptn.txt
    mkdir -p /media/glados/Data_Storage_2/hybrid/hybrid_mc/$dir/pass4/
    tlhbgeomnp_main.run -sd 1 -fd 6 -no_bw -f -i want_wsds.tlfptn.txt -o /media/glados/Data_Storage_2/hybrid/hybrid_mc/$dir/pass4/

    rm want_tamc_all.ok.txt
    rm want_wsds.txt

    echo "Running pass5 on $dir MC"
    mkdir -p /media/glados/Data_Storage_2/hybrid/hybrid_mc/$dir/pass5
    cd /media/glados/Data_Storage_2/hybrid/hybrid_mc/$dir/pass4/
    for file in `ls /media/glados/Data_Storage_2/hybrid/hybrid_mc/$dir/pass4/*`; do
	echo
	echo $file
	echo
	dev_stpfl12_main -det 37 -o /media/glados/Data_Storage_2/hybrid/hybrid_mc/$dir/pass5/${file:61:-27}.ps5.wsds.tlhbgeomnp.dst.gz -fit 7 $file  
    done
	
    
    cd $HYBRID_MC
done

echo "Concatinating DST files..."
dstcat.run */pass5/* -o hybrid_mc.all.dst.gz

echo "Pruning DST events that didn't trigger SD..."
dstsel.run +prfc,talex00 /home/glados/Data_Storage_2/hybrid/hybrid_mc/hybrid_mc.all.dst.gz -o /home/glados/Data_Storage_2/hybrid/hybrid_mc/hybrid_mc.all.cut.dst.gz

echo "Converting DST to ROOT file..."
dst2rt_ta.run /home/glados/Data_Storage_2/hybrid/hybrid_mc/hybrid_mc.all.cut.dst.gz 
