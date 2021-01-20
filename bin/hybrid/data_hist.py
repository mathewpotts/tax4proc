#!/usr/bin/env python3

# Import Libs
import ROOT
import TDSTio
import TLPY
import os
import sys
import argparse

# Call in external C script with quality cut functions
tax4proc_bin = os.environ['TAX4PROC_BIN']
ROOT.gROOT.LoadMacro("{0}/hybrid/qualct.C".format(tax4proc_bin))

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Concatinate data files (optional) and plot the data energy histogram.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    parser.add_argument('-concat',action='store_true',help='Use this flag if you want to concatinate all the data files into a single file.',default=False)
    parser.add_argument('-hyb',action='store_true',help='Use this flag if you want are dealing with hybrid data files.',default=False)
    args = parser.parse_args()
    concat=args.concat
    det=args.det
    hyb=args.hyb
    if det == 'mdtax4':
        fd      = 'TAX4_MIDDLE_DRUM'
        fd_path = os.environ['MDTAX4_DATA_PASS5']
    if det == 'brtax4':
        fd      = 'TAX4_BLACK_ROCK'
        fd_path = os.environ['BRTAX4_DATA_PASS5']
    return concat,fd,fd_path,det,hyb

def concat_dst(fd,fd_path,date,det):
    # Generate directories
    os.system("mkdir -p {0}/processing/concat_dst/".format(fd_path))
    os.system("mkdir -p {0}/processing/concat_drmdpcgf_dst/".format(fd_path))
    os.system("mkdir -p {0}/processing/root_all_files/".format(fd_path))
    
    # Concatinate files
    os.system("dstcat.run -o {0}/processing/concat_dst/20190625-{1}.{2}.ps5.ps2g.dst.gz {0}/[0-9]*[0-9]/*ps2g.dst.gz".format(fd_path,date,det))

    # Apply DCR to ps5 file
    os.system("dcrmdpcgf.run -f -fd {0} {1}/processing/concat_dst/20190625-{2}.{3}.ps5.ps2g.dst.gz -o {1}/processing/concat_drmdpcgf_dst/".format(fd,fd_path,date,det))
    
    # Convert the DST files into ROOT trees
    os.system("dst2rt_ta.run -f {0}/processing/concat_drmdpcgf_dst/20190625-{1}.{2}.ps5.ps2g.drmdpcgf.dst.gz -o {0}/processing/root_all_files/".format(fd_path,date,det))

def hybrid(hybrid_path,last_date,det_id):
    # Load pass5 root tree
    filepath = "{0}/pass5/hybrid_190625-{1}.tlhbgeomnp.root".format(fd_path,last_date[2:])
    file = ROOT.TFile(filepath,"r")
    taTree = file.Get("taTree")
    
    # Make TH1
    hEdata = ROOT.TH1D("hEdata",";log_{10}(E/eV);N_{EVENTS} / BIN",nxbins,log10en_lo,log10en_up);
    
    # Draw and save energy histogram
    taTree.Draw("log10((missing_E_corr(prfc.eng[{0}])))>>hEdata".format(det_id),"good_weather_cut(tlweat.code,fscn1->mir[1])==1 && log10(missing_E_corr(prfc.eng[{0}]))>=18.5".format(det_id))
    hEdata.SaveAs("{0}/hEdata.root".format(fd_path))

    ROOT.gPad.SaveAs("{0}/plots/hEdata_190625-{1}.png".format(fd_path,last_date[2:]))
    
    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')
    sys.exit(0)
    
    
def find_last_date(fd_path):
    file = os.popen('ls -drt {0}/[0-9][0-9]* | sort -n | tail -n1'.format(fd_path)).read()
    date = file[-9:-1]
    return date
    
if __name__=='__main__':
    # Define variables
    nxbins        = 20
    log10en_lo    = 18.0
    log10en_up    = 20.0
    concat,fd,fd_path,det,hyb = read_in_args()
    last_date = find_last_date(fd_path)
    if hyb is True:
            fd_path = os.environ['{0}_HYBRID_ROOT'.format(det.upper())]

    # Hybrid data?
    if hyb is True:
        if det == 'mdtax4':
            det_id = 6
        if det == 'brtax4':
            det_id = 7
        hybrid(fd_path,last_date,det_id)
    
    # Concatinate files?
    if concat is True:
        concat_dst(fd,fd_path,last_date,det)
    
    # Load pass5 root tree
    filepath = "{0}/processing/root_all_files/20190625-{1}.{2}.ps5.ps2g.drmdpcgf.root".format(fd_path,last_date,det)
    file = ROOT.TFile(filepath,"r")
    taTree = file.Get("taTree")
    
    # Make TH1
    hEdata = ROOT.TH1D("hEdata",";log_{10}(E/eV);N_{EVENTS} / BIN",nxbins,log10en_lo,log10en_up);
    
    # Draw and save energy histogram
    taTree.Draw("log10((missing_E_corr(prfc.eng[etrack.udata[0]])))>>hEdata","good_weather_cut(tlweat.code,fscn1->mir[1])==1 && etrack.qualct==1 && log10(missing_E_corr(prfc.eng[etrack.udata[0]]))>=18.5 && (prfc->chi2[etrack->udata[0]]/prfc->ndf[etrack->udata[0]])<20")
    hEdata.SaveAs("{0}/processing/hEdata.root".format(fd_path))

    plot_path = os.environ['{0}_DATA_ROOT'.format(det.upper())]
    ROOT.gPad.SaveAs("{0}/plots/hEdata_20190625-{1}.png".format(plot_path,last_date))
    
    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')
