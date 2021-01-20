#!/usr/bin/env python3

# Import Libs
import ROOT
import os
import argparse

# TAX4PROC DIR GLOBAL VARS
tax4proc_path     = os.environ['TAX4PROC']
tax4proc_bin_path = os.environ['TAX4PROC_BIN']

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Concatinate data files (optional) and plot the data energy histogram.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    args = parser.parse_args()
    det=args.det
    if det == 'mdtax4':
        fd_path = os.environ['MDTAX4_DATA_ROOT']
    if det == 'brtax4':
        fd_path = os.environ['BRTAX4_DATA_ROOT']
    return det,fd_path

def get_ontime(det,fd_path):
    filename    = '{0}/ontime_list_{1}.txt'.format(tax4proc_path,det)
    ontime_list = os.popen('readlink -f {0}/ontime/[0-9][0-9]*/*[0-9][0-9][0-9]*'.format(fd_path)).read().split('\n')
    del ontime_list[-1]
    for otf in ontime_list:
        if otf[-42:-18] == '20190625/y2019m06d25p001':
            #print(otf)
            os.system('readlink -f {0} > {1}'.format(otf,filename))
        elif int(otf[-42:-34]) >= 20190625:
            #print(otf)
            os.system('readlink -f {0} >> {1}'.format(otf,filename))
        else:
            continue
    os.system('{0}/calc_ontime.py -i {1} -det {2}'.format(tax4proc_bin_path,filename,det))
    with open(filename,'r') as f:
        lines = f.readlines()
    return float(lines[-1])

def find_last_date(fd_path):
    file = os.popen('ls -drt {0}/[0-9][0-9]* | sort -n | tail -n1'.format(fd_path)).read()
    date = file[-9:-1]
    return date

if __name__=='__main__':
    # Define variables
    det,fd_path = read_in_args()
    gwontime = get_ontime(det,fd_path)
    filepath = '{0}/mc/processing/hAp.root'.format(fd_path)

    # Open aperture root tree
    file = ROOT.TFile(filepath,"r")
    hAp = file.Get("hAp")

    # Scale histogram
    hExposure = hAp.Clone("hExposure")
    hExposure.Scale(gwontime * 3600)

    # Fit function
    func = ROOT.TF1("fit_func","TMath::Power(10,[0]*(1-TMath::Exp(([1]-x)/[2])))",18.,20.5)
    func.SetParameters(9,16,1)
    hExposure.Fit("fit_func","","",18.5,20.5)
    hExposure.GetYaxis().SetTitle("Exposure (m^{2}sr s)")
    hExposure.Draw()
    ROOT.gPad.SetLogy()
    ROOT.gPad.Update()
    hExposure.SaveAs("{0}/mc/processing/hExposure.root".format(fd_path))

    last_date = find_last_date(os.environ['{0}_DATA_PASS5'.format(det.upper())])
    ROOT.gPad.SaveAs("{0}/plots/hExposure_20190625-{1}.png".format(os.environ['{0}_DATA_ROOT'.format(det.upper())],last_date))
    
    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')
