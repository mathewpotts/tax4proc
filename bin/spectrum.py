#!/usr/bin/env python3

# Import Libs
import ROOT
import TDSTio
import os
import argparse
import numpy as np

tax4proc = os.environ['TAX4PROC']

def read_in_args():
    parser = argparse.ArgumentParser(description = 'Plot the aperture for detector.')
    parser.add_argument('-det',metavar='detector',action='store',help='Detector that you want to produce a data histogram for. (brtax4/mdtax4)',required=True)
    parser.add_argument('-hyb',action='store_true',help='Use this flag if you want are dealing with hybrid data files.',default=False)
    args = parser.parse_args()
    det  = args.det
    hyb  = args.hyb
    if det == 'mdtax4':
        fd_path = os.environ['MDTAX4_DATA_ROOT']
        hyb_path = os.environ['MDTAX4_HYBRID_ROOT']
    elif det == 'brtax4':
        fd_path = os.environ['BRTAX4_DATA_ROOT']
        hyb_path = os.environ['BRTAX4_HYBRID_ROOT']
    return det,fd_path,hyb,hyb_path

if __name__=='__main__':
    # Define variables
    det,fd_path,hyb,hyb_path = read_in_args()
    filepath_data            = '{0}/pass5/data/processing/hEdata.root'.format(fd_path) if not hyb else '{0}/hEdata.root'.format(hyb_path)
    filepath_exp             = '{0}/mc/processing/hExposure.root'.format(fd_path) if not hyb else '{0}/hExposure.root'.format(hyb_path)
    fldc_err                 = ROOT.TFeldmanCousins(0.65)
    
    # Grabbing Energy histogram for the data root file
    fEdata = ROOT.TFile(filepath_data,"r")
    hEdata = fEdata.Get("hEdata")

    # Grab variables from the data energy histogram
    nxbins = hEdata.GetNbinsX()
    log10en_lo = hEdata.GetXaxis().GetXmin()
    log10en_up = hEdata.GetXaxis().GetXmax()

    # Grabbing Exposure histogram and fitted function from the root file
    fExp = ROOT.TFile(filepath_exp,"r")
    hExposure = fExp.Get("hExposure")
    func_hExposure = hExposure.GetFunction("fit_func")

    # Defining TGraph
    gE3J = ROOT.TGraphAsymmErrors(1)
    gJ = ROOT.TGraphAsymmErrors(1)

    # Filling TGraph
    for i in range(1,nxbins):
        # Defining variables for each bin
        log_bin_center = hEdata.GetBinCenter(i)
        log_bin_width = hEdata.GetBinWidth(i)
        log_bin_lower_edge = (log_bin_center - log_bin_width / 2.0)
        log_bin_upper_edge = (log_bin_center + log_bin_width / 2.0)
        bin_center = pow(10,log_bin_center)
        bin_lower_edge = pow(10,log_bin_lower_edge)
        bin_upper_edge = pow(10,log_bin_upper_edge)

        if log_bin_center < 18.0:
            continue
        if log_bin_center > 21.5:
            continue

        # Each piece of the spectrum calculated for each bin
        nevent   = hEdata.GetBinContent(i)
        Ec3      = pow(bin_center,3)/1e24
        dE       = (bin_upper_edge - bin_lower_edge)
        exposure = (hExposure.GetBinContent(i+10))

        if nevent < 20.0 and nevent > 0.0:
            n_e1 = nevent - fldc_err.CalculateLowerLimit(nevent,0.0)
            n_e2 = fldc_err.CalculateUpperLimit(nevent,0.0) - nevent
        elif nevent >= 20:
            n_e1 = np.sqrt(nevent)
            n_e2 = n_e1
        else:
            n_e1 = 0
            n_e2 = 0

        gE3J.SetPoint(i-1, log_bin_center, Ec3 * nevent / dE / exposure)
        gE3J.SetPointError(i-1,0,0, Ec3 * n_e1 / dE / exposure, Ec3 * n_e2 / dE / exposure)

        gJ.SetPoint(i-1, log_bin_center, nevent / dE / exposure)
        gJ.SetPointError(i-1,0,0, n_e1 / dE / exposure, n_e2 / dE / exposure)   

    # Geting the TA combined spectrum data from root file
    ta_comp = ROOT.TFile("{0}/ta_combined_icrc2019.root".format(tax4proc),"r")
    compspec = ta_comp.Get("compspec")

    # Setting variables for combined spec TGraph 
    nentries = compspec.GetEntries()
    BC   = []
    FLUX = []
    E3J  = []
    kRed = 2
    
    # Fill arrays of values for the combined spec TGraph
    for i in range(0,nentries):
        compspec.GetEntry(i)
        BC.append(compspec.BinCenter)
        E = pow(10,compspec.BinCenter)
        J = compspec.J
        FLUX.append(J)
        E3 = pow(E,3)/1e24
        E3J.append(E3*J)

    flux = ROOT.TCanvas("flux","flux",900,600)
    flux.cd()
    flux.SetLogy()
    gJ.SetMarkerStyle(21)
    gJ.SetTitle(";log_{10}(E/eV);J(E) (eV^{-1}m^{-2}sr^{-1}s^{-1})")
    gJ.GetYaxis().SetRangeUser(1e-34,2e-30)
    gJ.Draw("a,e1p")

    icrc_compspec = ROOT.TGraph(nentries,np.array(BC),np.array(FLUX))
    icrc_compspec.SetMarkerColor(kRed)
    icrc_compspec.SetMarkerStyle(20)
    icrc_compspec.Draw("P")

    # Setting up legend
    legend = ROOT.TLegend(0.1,0.7,0.3,0.9)
    legend.SetHeader("{0} Monocular Spectrum".format(det.upper()),"C")
    icrcentry = legend.AddEntry("icrc_flux","ICRC2019 Combined Spectrum","p")
    icrcentry.SetMarkerColor(2)
    icrcentry.SetMarkerStyle(20)
    dataentry = legend.AddEntry("gJ","{0} Monocular Spectrum".format(det.upper()),"p")
    dataentry.SetMarkerStyle(21)
    legend.Draw()
    
    ROOT.gPad.Update()
    print('\nPress enter to continue.\n')
    os.system('read gonext')
    
    spec = ROOT.TCanvas("spec","spec",900,600)
    spec.cd()
    spec.SetLogy()
    gE3J.SetMarkerStyle(21)
    gE3J.SetTitle(";log_{10}(E/eV);E^{3}*J(E) / 1x10^{24} (eV^{2}m^{-2}sr^{-1}s^{-1})")
    gE3J.GetYaxis().SetRangeUser(0.1,10.0)
    gE3J.Draw("a,e1p")

    # Setting up TGraph
    icrc_compspec = ROOT.TGraph(nentries,np.array(BC),np.array(E3J))
    icrc_compspec.SetMarkerColor(kRed)
    icrc_compspec.SetMarkerStyle(20)
    icrc_compspec.Draw("P")

    # Setting up legend
    legend2 = ROOT.TLegend(0.1,0.7,0.3,0.9)
    legend2.SetHeader("{0} Monocular Spectrum".format(det.upper()),"C")
    icrcentry2 = legend2.AddEntry("icrc_compspec","ICRC2019 Combined Spectrum","p")
    icrcentry2.SetMarkerColor(2)
    icrcentry2.SetMarkerStyle(20)
    dataentry2 = legend2.AddEntry("gE3J","{0} Monocular Spectrum".format(det.upper()),"p")
    dataentry2.SetMarkerStyle(21)
    legend2.Draw()

    # Save spectrum and flux
    ROOT.gPad.Update()
    if not hyb:
        flux.SaveAs('{0}/plots/monocular_flux.png'.format(fd_path))
        spec.SaveAs('{0}/plots/monocular_spectrum.png'.format(fd_path))
    else:
        flux.SaveAs('{0}/plots/monocular_flux.png'.format(hyb_path))
        spec.SaveAs('{0}/plots/monocular_spectrum.png'.format(hyb_path))
    
    # Prompt exit
    print('\nPress enter to continue.\n')
    os.system('read gonext')

