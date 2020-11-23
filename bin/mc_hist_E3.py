#!/usr/bin/env python3

# Import Libs
import ROOT
import TDSTio
import TLPY
import os
import argparse

# Call in external C script with quality cut functions
ROOT.gROOT.LoadMacro("hybrid/qualct.C")

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Concatinate data files (optional) and plot the mc energy histogram.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    parser.add_argument('-concat',action='store_true',help='Use this flag if you want to concatinate all the data files into a single file.',default=False)
    args = parser.parse_args()
    concat=args.concat
    det=args.det
    if det == 'mdtax4':
        fd = 'TAX4_MIDDLE_DRUM'
        fd_path = os.environ['MDTAX4_MC_ROOT']
    if det == 'brtax4':
        fd = 'TAX4_BLACK_ROCK'
        fd_path = os.environ['BRTAX4_MC_ROOT']
    return concat,fd,fd_path,det

def concat_dst(fd,fd_path,date,det):
    # Generate directories
    os.system("mkdir -p {0}/processing/concat_dst/".format(fd_path))
    os.system("mkdir -p {0}/processing/concat_drmdpcgf_dst/".format(fd_path))
    os.system("mkdir -p {0}/processing/root_all_files/".format(fd_path))
        
    # Concatinate files
    os.system("dstcat.run -o {0}/processing/concat_dst/pr.1e17-1e21.20190625-{1}.{2}.ps5.dst.gz {0}/*/pass5/*.dst".format(fd_path,date,det))
    os.system("dstcat.run -o {0}/processing/concat_dst/pr.1e17-1e21.20190625-{1}.{2}.ps0.dst.gz {0}/*/pass0/*.dst".format(fd_path,date,det))

    # Apply DCR to ps5 file
    os.system("dcrmdpcgf.run -f -fd {0} {1}/processing/concat_dst/pr.1e17-1e21.20190625-{2}.{3}.ps5.dst.gz -o {1}/processing/concat_drmdpcgf_dst/".format(fd,fd_path,date,det))
    # Convert the DST files into ROOT trees
    os.system("dst2rt_ta.run -f {0}/processing/concat_drmdpcgf_dst/pr.1e17-1e21.20190625-{1}.{2}.ps5.drmdpcgf.dst.gz -o {0}/processing/root_all_files/".format(fd_path,date,det))
    os.system("dst2rt_ta.run -f {0}/processing/concat_dst/pr.1e17-1e21.20190625-{1}.{2}.ps0.dst.gz -o {0}/processing/root_all_files/".format(fd_path,date,det))

def find_last_date(fd_path):
    file = os.popen('ls -drt {0}/[0-9][0-9]* | sort -n | tail -n1'.format(fd_path)).read()
    date = file[-9:-1]
    return date


if __name__=='__main__':
    # Define variables
    concat,fd,fd_path,det = read_in_args()
    last_date = find_last_date(fd_path)
    
    # Concatinate files?
    if concat is True:
        concat_dst(fd,fd_path,last_date,det)

    # Create Canvas
    canvas = ROOT.TCanvas("c1")
    canvas.cd()
    
    # Load mc pass5 root tree
    filepath = "{0}/processing/root_all_files/pr.1e17-1e21.20190625-{1}.{2}.ps5.drmdpcgf.root".format(fd_path,last_date,det)
    print("Opening..", filepath)
    file = ROOT.TFile(filepath,"r")
    taTree = file.Get("taTree")
    
    # Make TH1
    hErecon  = ROOT.TH1D("hErecon",";log_{10}(E/eV);N_{EVENTS} / BIN",40,17,21)
    
    
    # Draw Recon and save energy histogram
    taTree.Draw("log10((missing_E_corr(prfc->eng[etrack->udata[0]])))>>hErecon","(etrack->qualct==1) * weight(log10(mc04->energy))","sames")
    hErecon.SaveAs("{0}/processing/hErecon_E3.root".format(fd_path))
    ROOT.gPad.Update()
    
    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')
    print('Plotting thrown...')

    # Load mc pass0 root tree
    filepath = "{0}/processing/root_all_files/pr.1e17-1e21.20190625-{1}.{2}.ps0.root".format(fd_path,last_date,det)
    file = ROOT.TFile(filepath,"r")
    taTree2 = file.Get("taTree")

    # Draw thrown and save energy histogram
    hEthrown = ROOT.TH1D("hEthrown",";log_{10}(E/eV);N_{EVENTS} / BIN",40,17,21)
    taTree2.Draw("log10(mc04->energy)>>hEthrown", "weight(log10(mc04->energy))")
    hEthrown.SaveAs("{0}/processing/hEthrown_E3.root".format(fd_path))
    ROOT.gPad.Update()
    
    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')
