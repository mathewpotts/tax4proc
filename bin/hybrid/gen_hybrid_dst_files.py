#!/usr/bin/env python3

# Import Libs
import TLPY
import TDSTio
import argparse
from multiprocessing import Pool,cpu_count
import os
import time

starttime = time.time()

# Grab environmental variables from bash
sdtax4     = os.environ['SDTAX4_DATA_ROOT']
listfiles  = sdtax4+'/utsdpass0/listfiles'
dstfiles   = sdtax4+'/utsdpass0/dst'


def read_in_args():
    parser = argparse.ArgumentParser(description = 'Concatinate data files (optional) and plot the data energy histogram.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    parser.add_argument('-prep',action='store_true',help='Use this flag if you want to concatinate all the data files into a single file.',default=False)
    args = parser.parse_args()
    prep=args.prep
    det=args.det
    if det == 'mdtax4':
        hyb_path   = os.environ['MDTAX4_HYBRID_ROOT']
        fd_path    = os.environ['MDTAX4_DATA_ROOT']
    if det == 'brtax4':
        hyb_path   = os.environ['BRTAX4_HYBRID_ROOT']
        fd_path    = os.envrion['BRTAX4_DATA_ROOT']
    return prep,hyb_path,fd_path,det 

class gen_hybrid_dst_files:
    def __init__(self,START,END,det):
        os.system('''
        mkdir -p {1}/mc
        mkdir -p {1}/tmatch
        mkdir -p {1}/fd_byday
        mkdir -p {0}/utsdpass0/listfiles
        mkdir -p {0}/utsdpass0/dst
        '''.format(sdtax4,hyb_path))
        self.start    = START
        self.end      = END
        self.hyb_path = hyb_path 
        
    def prep_listfile_parallel(self,tower_path):
        os.system('utsdpass0_prep_listfiles.py -d1 {1} -d2 {2} -i {0}/DATA -o {0}'.format(tower_path,self.start,self.end))

    def prep_cleanup(self,tower):
        for i in range(self.start,self.end):
            dstfile="want_tasd_{0}.ok.txt".format(i)
            if os.path.isfile(tower+'/'+dstfile):
                want_file = tower+'/'+dstfile
                os.system('cat {0} >> {1}/{2}'.format(want_file,listfiles,dstfile))
        os.system('find {0} -empty -type f -delete'.format(listfiles))
        
    def utsdpass0(self,listfile):
        date = listfile[-14:-8]
        print('Processing {0}...'.format(listfile[:-1]))
        os.system('utsdpass0.run -f -incomplete -o {1} -d {2} -i {0} > {1}/SD{2}.log'.format(listfile[:-1],dstfiles,date))
        
    def fd_byday(self,fd_path):
        fd_byday = self.hyb_path+'/fd_byday'
        os.system('readlink -f {0}/pass3/data/*/*ps2g.dst > {1}/want_tafd_all.ok.txt'.format(fd_path,fd_byday))
        os.system('ta_byday.run -i {0}/want_tafd_all.ok.txt -f -o {0} -v'.format(fd_byday))
            
    def tmatch(self,fd_day):
        date       = fd_day[-14:-8]
        tmatch_dir = self.hyb_path+'/tmatch' 
        if os.path.isfile('{0}/SD{1}.incomplete.sdpass0.dst.gz'.format(dstfiles,date)):
            print('tmatch.run -ist 2 -f -matched_only -ifl2 {1}/SD{2}.incomplete.sdpass0.dst.gz -ifl1 {0}'.format(fd_day,dstfiles,date))
            os.system('tmatch.run -ist 2 -f -matched_only -ifl2 {1}/SD{2}.incomplete.sdpass0.dst.gz -ifl1 {0}'.format(fd_day,dstfiles,date))
            os.system('mv tmatch_{0}*.gz {1}'.format(date,tmatch_dir))
            
        
#    def pass5_recon(self):
#        dev_stpfl12_main -det 37 -db -fit 7 -xm 685,720,755,790,825,860 $TMATCH/tmatch_mo_int2_all.tlhbgeomnp.dst.gz -o $TMATCH/pass5/hybrid_${START}_${END}.tlhbgeomnp.dst.gz > $TMATCH/pass5/hybrid_${START}_${END}.tlhbgeomnp.log
#        dst2rt_ta.run -f $TMATCH/pass5/hybrid_${START}_${END}.tlhbgeomnp.dst.gz -o $TMATCH/pass5
        
if __name__ == '__main__':
    prep,hyb_path,fd_path,det = read_in_args()
    gen = gen_hybrid_dst_files(190625,200924,det)

    # Raw SD directories
    bf = sdtax4+"/raw/bf"
    dm = sdtax4+"/raw/dm"
    km = sdtax4+"/raw/km"
    sc = sdtax4+"/raw/sc"
    sn = sdtax4+"/raw/sn"
    sr = sdtax4+"/raw/sr"

    # Pick towers in proper lobe
    if det == 'mdtax4':
        towers = [dm,km,sn]
    elif det == 'brtax4':
        towers = [bf,sr,sc]
    else:
        exit(1)

    if prep:
        # Run 6 instances of the list prep function
        with Pool(cpu_count()) as p:
            #p.map(gen.prep_listfile_parallel,towers)
            
            # Move and clean up prep files
            #p.map(gen.prep_cleanup,towers)
            
            files = list(os.popen('ls {0}/*'.format(listfiles)))
            p.map(gen.utsdpass0,files)
        os.system('mv *.dst.gz {0}'.format(dstfiles))

    #gen.fd_byday(fd_path)

    #fd_day_list = list(os.popen('ls {0}/fd_byday/*dst.gz'.format(hyb_path)))
    #print(fd_day_list[1])
    #with Pool(cpu_count()) as p:
    #    p.map(gen.tmatch,fd_day_list)
    

print("-- %s seconds ---" % (time.time() - starttime ))