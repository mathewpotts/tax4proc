#!/usr/bin/env python2

# Import Libs
import ROOT
import os
import argparse
import numpy as np

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Plot the aperture for detector.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    args = parser.parse_args()
    det=args.det
    return det

if __name__=='__main__':
    # Define variables
    det             = read_in_args()
    filepath_recon  = '/home/glados/Data_Storage_2/hybrid/{0}/hybrid_mc/processing/hErecon_E3.root'.format(det)
    filepath_thrown = '/home/glados/Data_Storage_2/{0}/mc/processing/hEthrown_E3.root'.format(det)
    rpmin           = 100.0
    rpmax           = 50.0e3
    theta_max       = 70.0
    area            = np.pi * (rpmax * rpmax - rpmin * rpmin)
    omega           = 2 * np.pi * (1 - np.cos(np.deg2rad(theta_max)))
    area_x_omega    = area * omega
    
    
    # Grab recon and thrown root trees
    file1 = ROOT.TFile(filepath_recon,"r")
    hErecon = file1.Get("hErecon")
    hErecon.Sumw2()
    file2 = ROOT.TFile(filepath_thrown,"r")
    hEthrown = file2.Get("hEthrown")
    hEthrown.Sumw2()

    # Make hAp plot
    hAp = hErecon.Clone("hAp")
    hAp.Sumw2()
    hAp.Divide(hEthrown)
    hAp.Scale(area_x_omega)
    hAp.GetYaxis().SetTitle("Reconstructed Aperture (m^{2}sr)")
    
    # Fit a function to aperture plot
    func = ROOT.TF1("fit_func","TMath::Power(10,[0]*(1-TMath::Exp(([1]-x)/[2])))",18.,20.5)
    func.SetParameters(9,16,1)
    hAp.Fit("fit_func","","",18.5,20.5)

    # Draw hAp
    hAp.Draw()
    ROOT.gPad.SetLogy()
    hAp.SaveAs("/home/glados/Data_Storage_2/{0}/mc/processing/hAp.root".format(det))
    
    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')
