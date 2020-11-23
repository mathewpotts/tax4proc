#!/usr/bin/env python2

# Import Libs
import ROOT
import os
import argparse

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Concatinate data files (optional) and plot the data energy histogram.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    args = parser.parse_args()
    det=args.det
    return det

def get_ontime(det):
    filename = '/home/glados/Data_Storage_2/scripts/ontime_list_{0}.txt'.format(det)
    os.system('readlink -f /media/glados/Data_Storage_2/{0}/ontime/*/* > {1}'.format(det,filename))
    os.system('/home/glados/Data_Storage_2/scripts/calc_ontime.py -i {0} -det {1}'.format(filename,det))
    with open(filename,'r') as f:
        lines = f.readlines()
    return float(lines[-1])

if __name__=='__main__':
    # Define variables
    det = read_in_args()
    gwontime = get_ontime(det)
    filepath = '/media/glados/Data_Storage_2/{0}/mc/processing/hAp.root'.format(det)

    # Open aperture root tree
    file = ROOT.TFile(filepath,"r")
    hAp = file.Get("hAp")

    # Scale histogram
    hExposure = hAp.Clone("hExposure")
    hExposure.Scale(gwontime * 3600)

    # Fit function
    func = ROOT.TF1("fit_func","TMath::Power(10,[0]*(1-TMath::Exp(([1]-x)/[2])))",18.,20.5)
    func.SetParameters(9,16,1)
    hExposure.Fit("fit_func","","",18.0,20.5)
    hExposure.GetYaxis().SetTitle("Exposure (m^{2}sr s)")
    hExposure.Draw()
    hExposure.SaveAs("/media/glados/Data_Storage_2/{0}/mc/processing/hExposure.root".format(det))

    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')
