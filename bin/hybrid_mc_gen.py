#!/usr/bin/env python3
import os
import sys
import argparse
from multiprocessing import Pool

class hybrid_mc_gen:
    def __init__(self,DET):
        if det == 'mdtax4':
            DET_ID = 37
            self.FD_ID  = 6
            self.SD_ID  = 1
        elif det == 'brtax4':
            DET_ID = 38
            self.FD_ID  = 7
            self.SD_ID  = 2
        else:
            print("{0} isn't mdtax4 or brtax4. Try again.".format(DET))
            sys.exit(1)
        self.HYBRID_MC_PATH = os.environ['{0}_HYBRID_ROOT'.format(DET.upper())] + '/mc'
        self.MONO_MC_PATH   = os.environ['{0}_MC_ROOT'.format(DET.upper())]
        print('''
        Generating hybrid MC for {0}.
        Setting ID(s)
        DET_ID : {1}
        FD_ID  : {2}
        SD_ID  : {3}
        '''.format(DET.upper(),DET_ID,self.FD_ID,self.SD_ID))
    
    def mkdirs(self,DATE):
        os.system('''
        mkdir -p {0}/{1}/pass3/tl_1sd_hy_sim
        mkdir -p {0}/{1}/pass4/
        mkdir -p {0}/{1}/pass5/
        '''.format(self.HYBRID_MC_PATH,DATE))

    def tl_1sd_hy_simulation(self,DATE):
        print("Adding SD simulation data to {0} FD MC...".format(DATE))
        os.system('''
        readlink -f {1}/{2}/pass3/*ps3.dst > {0}/{2}/want_tamc_all.ok.txt; 
        tl_1sd_hy_simulation.run -f -sd_det {3} -i {0}/{2}/want_tamc_all.ok.txt; 
        mv *wsds.dst.gz {0}/{2}/pass3/tl_1sd_hy_sim/; 
        rm -rvf {0}/{2}/want_tamc_all.ok.txt;
        '''.format(self.HYBRID_MC_PATH,self.MONO_MC_PATH,DATE,self.SD_ID))

    def tlfptn(self,DATE):
        print("Adding tlfptn bank to {0} MC".format(DATE))
        os.system('''
        readlink -f {0}/{1}/pass3/tl_1sd_hy_sim/*ps3.wsds.dst.gz > {0}/{1}/want_wsds.txt;
        tlfptn.run -det {2} -i {0}/{1}/want_wsds.txt -f -no_bw -o {0}/{1}/pass3;
        rm -rvf {0}/{1}/want_wsds.txt;
        '''.format(self.HYBRID_MC_PATH,DATE,self.SD_ID))

    def tlhbgeonp(self,DATE):
        print("Running tlhbgeonp on {0} MC".format(DATE))
        os.system('''
        readlink -f {0}/{1}/pass3/*.tlfptn.dst.gz > {0}/{1}/want_wsds.tlfptn.txt;
        tlhbgeomnp_main.run -sd {3} -fd {2} -no_bw -f -i {0}/{1}/want_wsds.tlfptn.txt -o {0}/{1}/pass4;
        rm -rvf {0}/{1}/want_wsds.tlfptn.txt
        '''.format(self.HYBRID_MC_PATH,DATE,self.FD_ID,self.SD_ID))

    def multi_pass5(self,DATE):
        print("Running pass5 on {0} MC".format(DATE))
        PARTS = os.popen('readlink -f {0}/{1}/pass4/*'.format(self.HYBRID_MC_PATH,DATE)).read().split('\n')
        del PARTS[-1] # deleting last list entry -> ''
        with Pool(24) as p:
            p.map(dev_stpfl12_main,PARTS)
            pass

    def concat(self):
        os.system('''
        echo "Concatinating DST files...";
        dstcat.run {0}/*/pass5/* -o {0}/hybrid_mc.all.dst.gz;
        echo "Pruning DST events that didn't trigger SD...";
        dstsel.run -f +prfc,talex00 {0}/hybrid_mc.all.dst.gz -o {0}/hybrid_mc.all.cut.dst.gz;
        echo "Converting DST to ROOT file...";
        dst2rt_ta.run -f {0}/hybrid_mc.all.cut.dst.gz;
        '''.format(self.HYBRID_MC_PATH))

def dev_stpfl12_main(PART):
    DET          = PART[-80:-74]
    DET_ID       = 37 if DET == 'mdtax4' else 38 
    file_prefix  = PART[-55:-27]
    outfile_path = PART[:-61] + '/pass5/'
    os.system('''
    dev_stpfl12_main -det {3} -db -o {2}/{1}.ps5.wsds.tlhbgeomnp.dst.gz -fit 7 {0};
    '''.format(PART,file_prefix,outfile_path,DET_ID))

def main(DATES,DET):
    hyb_gen = hybrid_mc_gen(DET)
    for f in DATES:
        DATE = f[-8:]
        hyb_gen.mkdirs(DATE)
        hyb_gen.tl_1sd_hy_simulation(DATE)
        hyb_gen.tlfptn(DATE)
        hyb_gen.tlhbgeonp(DATE)
        hyb_gen.multi_pass5(DATE)
    hyb_gen.concat()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Concatinate data files (optional) and plot the data energy histogram.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a hybrid MC for. (brtax4/mdtax4)',required=True)
    args = parser.parse_args()
    det  = args.det
    MONO_MC_PATH   = os.environ['{0}_MC_ROOT'.format(det.upper())]
    DATES = os.popen('ls -d {0}/[0-9][0-9]*'.format(MONO_MC_PATH)).read().split('\n')
    del DATES[-1] # deleting last list entry -> ''
    main(DATES,det)

