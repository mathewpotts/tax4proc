#!/usr/bin/env python3

# Import Libs
import TLPY
import TDSTio
import argparse
from multiprocessing import Pool,cpu_count
import os
import subprocess
import time


starttime = time.time()

# Grab environmental variables from bash
sdtax4     = os.environ['SDTAX4_DATA_ROOT']
listfiles  = sdtax4+'/utsdpass0/listfiles'
dstfiles   = sdtax4+'/utsdpass0/dst'


def read_in_args():
    parser = argparse.ArgumentParser(description = "Used to generate hybrid data files by time matching events with Dmitri's tmatch program.")
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce hybrid data files for. (brtax4/mdtax4)',required=True)
    parser.add_argument('-prep',action='store_true',help='This flag preps the SD files for utsdpass0 in a date range',default=False)
    args = parser.parse_args()
    prep=args.prep
    det=args.det
    if det == 'mdtax4':
        hyb_path   = os.environ['MDTAX4_HYBRID_ROOT']
        fd_path    = os.environ['MDTAX4_DATA_ROOT']
        det_id     = 'TAX4_MIDDLE_DRUM'
        det_num    = 37
        fit    = 7
    if det == 'brtax4':
        hyb_path   = os.environ['BRTAX4_HYBRID_ROOT']
        fd_path    = os.envrion['BRTAX4_DATA_ROOT']
        det_id     = 'TAX4_BLACK_ROCK'
        det_num    = 38
        fit    = 8
    return prep,hyb_path,fd_path,det,det_id,det_num,fit 

class gen_hybrid_dst_files:
    def __init__(self,START,END,det):
        os.system('''
        mkdir -p {1}/mc
        mkdir -p {1}/tmatch
        mkdir -p {1}/fd_byday
        mkdir -p {1}/pass5
        mkdir -p {0}/utsdpass0/listfiles
        mkdir -p {0}/utsdpass0/dst
        '''.format(sdtax4,hyb_path))
        self.start    = START
        self.end      = END
        self.hyb_path = hyb_path
        
    def prep_listfile_parallel(self,tower_path):
        cmd = 'utsdpass0_prep_listfiles.py -d1 {1} -d2 {2} -i {0}/DATA -o {0}'.format(tower_path,self.start,self.end)
        print(cmd)
        subprocess.run(cmd.split(' '))

    def prep_cleanup(self,tower):
        for i in range(self.start,self.end):
            dstfile="want_tasd_{0}.ok.txt".format(i)
            if os.path.isfile(tower+'/'+dstfile):
                want_file = tower+'/'+dstfile
                os.system('cat {0} >> {1}/{2}'.format(want_file,listfiles,dstfile))
        os.system('find {0} -empty -type f -delete'.format(listfiles))
        
    def utsdpass0(self,listfile):
        date = listfile[-14:-8]
        cmd1 = 'utsdpass0.run -f -incomplete -o {1} -d {2} -i {0} > {1}/SD{2}.log'.format(listfile[:-1],dstfiles,date)
        cmd2 = 'tlfptn.run -f {0}/SD{1}.incomplete.sdpass0.dst.gz -o {0}'.format(dstfiles,date)
        print(cmd1)
        subprocess.run(cmd1.split(' '))
        print(cmd2)
        subprocess.run(cmd2.split(' '))
        
    # NEED TO FIX THIS FUNCTION TO WORK WITHOUT WANT FILE OR WITH CUSTOM WANT FILE    
    def fd_byday(self,fd_path):
        fd_byday = self.hyb_path+'/fd_byday'
        cmd1     = 'readlink -f {0}/pass3/data/*/*ps2g.dst* > {1}/want_tafd_all.ok.txt'.format(fd_path,fd_byday)
        cmd2     =  'ta_byday.run -i {0}/want_tafd_all.ok.txt -f -o {0} -v'.format(fd_byday) # only works if the want list is small
        print(cmd1)
        os.system(cmd1)
        print(cmd2)
        os.system(cmd2)
            
    def tmatch(self,fd_day):
        date       = fd_day[-14:-8]
        tmatch_dir = self.hyb_path+'/tmatch' 
        if os.path.isfile('{0}/SD{1}.incomplete.tlfptn.dst.gz'.format(dstfiles,date)):
            cmd1 = 'tmatch.run -f -ist 2 -matched_only -ifl1 {0} -ifl2 {1}/SD{2}.incomplete.tlfptn.dst.gz'.format(fd_day[:-1],dstfiles,date) 
            print(cmd1)
            os.system(cmd1)
            if os.path.isfile('tmatch_{1}.dst.gz'.format(tmatch_dir,date)):
                os.system('mv tmatch_{0}*.gz {1}'.format(date,tmatch_dir))
                cmd2 = 'tlfptn.run -f {0}/tmatch_{1}.dst.gz -o {0}'.format(tmatch_dir,date)
                print(cmd2)
                os.system(cmd2)
        
    def pass5_recon(self,det_num,fit):
        pass5_dir  = self.hyb_path + '/pass5'
        tmatch_dir = self.hyb_path+'/tmatch'
        cmd1       = 'dev_stpfl12_main -det {4} -db -fit {5} -xm 685,720,755,790,825,860 {3}/tmatch_mo_int2_all.tlhbgeomnp.dst.gz -o {2}/hybrid_{0}_{1}.tlhbgeomnp.dst.gz > {2}/hybrid_{0}_{1}.tlhbgeomnp.log'.format(self.start,self.end,pass5_dir,tmatch_dir,det_num,fit)
        cmd2       = 'dst2rt_ta.run -f {2}/hybrid_{0}_{1}.tlhbgeomnp.dst.gz -o {2}'.format(self.start,self.end,pass5_dir,tmatch_dir)
        print(cmd1)
        subprocess.run(cmd1.split(' '))
        print(cmd2)
        subprocess.run(cmd2.split(' '))

def cmd_run(cmd):
    print(cmd)
    subprocess.run(cmd.split(' '))
    
        
if __name__ == '__main__':
    prep,hyb_path,fd_path,det,det_id,det_num,fit = read_in_args()
    gen = gen_hybrid_dst_files(190625,201207,det)

    # Raw SD directories
    bf = sdtax4+"/raw/bf"
    dm = sdtax4+"/raw/dm"
    km = sdtax4+"/raw/km"
    sc = sdtax4+"/raw/sc"
    sn = sdtax4+"/raw/sn"
    sr = sdtax4+"/raw/sr"

    # Pick towers in proper lobe
    if det == 'mdtax4':
        towers = [dm,km,sn,bf,sr,sc]
    elif det == 'brtax4':
        towers = [dm,km,sn,bf,sr,sc]
    else:
        exit(1)

    if prep:
        # Run 8 instances of the list prep function
        
        with Pool(cpu_count()) as p:
            #p.map(gen.prep_listfile_parallel,towers)             # Generate listfiles of SD files by date
            #p.map(gen.prep_cleanup,towers)                       # Move and clean up prep files
            files = list(os.popen('ls {0}/*'.format(listfiles))) # Create a list of listfiles for next step
            p.map(gen.utsdpass0,files)                           # Run utsdpass0 on each listfile
        #os.system('rm *.Y')

        #gen.fd_byday(fd_path) # This script has a problem with large date ranges... need to fix

    fd_day_list = list(os.popen('ls {0}/fd_byday/*dst.gz'.format(hyb_path)))
    #with Pool(1) as p:
    #    p.map(gen.tmatch,fd_day_list)

    #concat all tmatch files
    cmd = 'dstcat.run -f {0}/tmatch/*[0-9].tlfptn.dst.gz -o {0}/tmatch/tmatch_mo_int2_all.tlfptn.dst.gz'.format(hyb_path)
    #print(cmd)
    #os.system(cmd)

    # pass 4
    #cmd_run('tlhbgeomnp_main.run -f -sd TASD_ALL -fd {0} {1}/tmatch/tmatch_mo_int2_all.tlfptn.dst.gz'.format(det_id,hyb_path))

    #cmd = 'mv tmatch_mo_int2_all.tlhbgeomnp.dst.gz {1}/tmatch/'.format(det_id,hyb_path)
    #print(cmd)
    #subprocess.run(cmd.split(' '))

    # run pass5
    #gen.pass5_recon(det_num,fit)

    
    

print("-- %s seconds ---" % (time.time() - starttime ))
