#!/usr/bin/env python3

# Import Libs
import ROOT
import os
import argparse
import numpy as np

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Plot the aperture for detector.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    parser.add_argument('-hyb',action='store_true',help='Use this flag if you want are dealing with hybrid data files.',default=False)
    args = parser.parse_args()
    det  = args.det
    hyb  = args.hyb
    if det == 'mdtax4':
        fd_path  = os.environ['MDTAX4_MC_ROOT']
        hyb_path = os.environ['MDTAX4_HYBRID_ROOT']
    if det == 'brtax4':
        fd_path = os.environ['BRTAX4_MC_ROOT']
        hyb_path = os.environ['BRTAX4_HYBRID_ROOT']
    return det,fd_path,hyb,hyb_path

def find_last_date(fd_path):
    file = os.popen('ls -drt {0}/[0-9][0-9]* | sort -n | tail -n1'.format(fd_path)).read()
    date = file[-9:-1]
    return date

if __name__=='__main__':
    # Define variables
    det,fd_path,hyb,hyb_path = read_in_args()
    filepath_recon           = '{0}/processing/hErecon_E3.root'.format(fd_path) if not hyb else '{0}/hErecon_E3.root'.format(hyb_path)
    filepath_thrown          = '{0}/processing/hEthrown_E3.root'.format(fd_path) if not hyb else '{0}/hEthrown_E3.root'.format(hyb_path)
    print(filepath_recon,filepath_thrown)
    rpmin                    = 100.0
    rpmax                    = 50.0e3
    theta_max                = 70.0
    area                     = np.pi * (rpmax * rpmax - rpmin * rpmin)
    omega                    = 2 * np.pi * (1 - np.cos(np.deg2rad(theta_max)))
    area_x_omega             = area * omega
    
    
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
    func = ROOT.TF1("fit_func","TMath::Power(10,[0]*(1-TMath::Exp(([1]-x)/[2])))",18.5,20.5)
    func.SetParameters(9,16,1)
    hAp.Fit("fit_func","","",18.5,20.5)

    # Draw hAp
    hAp.Draw()
    ROOT.gPad.SetLogy()
    ROOT.gPad.Update()
    if not hyb:
        hAp.SaveAs("{0}/processing/hAp.root".format(fd_path))
        plot_path = os.environ['{0}_DATA_ROOT'.format(det.upper())]
        ROOT.gPad.SaveAs("{0}/plots/hAperture_20190625-{1}.png".format(plot_path,find_last_date(os.environ['{0}_DATA_PASS5'.format(det.upper())])))
    else:
        hAp.SaveAs("{0}/hAp.root".format(hyb_path))
        ROOT.gPad.SaveAs("{0}/plots/hAperture_20190625-{1}.png".format(hyb_path,find_last_date(os.environ['{0}_DATA_PASS5'.format(det.upper())])))
    
    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')
