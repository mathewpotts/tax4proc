#!/usr/bin/env python3

# Import Lib
import ROOT
import TLPY
import TDSTio
import os
import numpy as np

# Define variables
det             = ['mdtax4','brtax4']
filepath_data_md   = '{0}/processing/hEdata.root'.format(os.environ['{0}_DATA_PASS5'.format(det[0].upper())])
filepath_data_br   = '{0}/processing/hEdata.root'.format(os.environ['{0}_DATA_PASS5'.format(det[1].upper())])
filepath_exp_md    = '{0}/processing/hExposure.root'.format(os.environ['{0}_MC_ROOT'.format(det[0].upper())])
filepath_exp_br    = '{0}/processing/hExposure.root'.format(os.environ['{0}_MC_ROOT'.format(det[1].upper())])
fldc_err = ROOT.TFeldmanCousins(0.65)
    
# Grabbing Energy histogram for the data root file
fEdata_md = ROOT.TFile(filepath_data_md,"r")
fEdata_br = ROOT.TFile(filepath_data_br,"r")
hEdata_md = fEdata_md.Get("hEdata")
hEdata_br = fEdata_br.Get("hEdata")

# Grab variables from the data energy histogram
nxbins_md = hEdata_md.GetNbinsX()
nxbins_br = hEdata_br.GetNbinsX()
log10en_lo_md = hEdata_md.GetXaxis().GetXmin()
log10en_lo_br = hEdata_br.GetXaxis().GetXmin()
log10en_up_md = hEdata_md.GetXaxis().GetXmax()
log10en_up_br = hEdata_br.GetXaxis().GetXmax()

# Grabbing Exposure histogram and fitted function from the root file
fExp_md = ROOT.TFile(filepath_exp_md,"r")
fExp_br = ROOT.TFile(filepath_exp_md,"r")
hExposure_md = fExp_md.Get("hExposure")
hExposure_br = fExp_br.Get("hExposure")
func_hExposure_md = hExposure_md.GetFunction("fit_func")
func_hExposure_br = hExposure_br.GetFunction("fit_func")

# Defining TGraph
gE3J = ROOT.TGraphAsymmErrors(20)
gJ = ROOT.TGraphAsymmErrors(20)

# Filling TGraph
for i in range(1,nxbins_md):
    # Defining variables for each bin
    log_bin_center = hEdata_md.GetBinCenter(i)
    log_bin_width = hEdata_md.GetBinWidth(i)
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
    nevent_md = hEdata_md.GetBinContent(i)
    nevent_br = hEdata_br.GetBinContent(i)
    nevent    = nevent_md + nevent_br
    #print(nevent)
    Ec3 = pow(bin_center,3)/1e24
    dE = (bin_upper_edge - bin_lower_edge)
    exposure_md = (hExposure_md.GetBinContent(i+10))
    exposure_br = (hExposure_br.GetBinContent(i+10))
    exposure    = exposure_md + exposure_br

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
    gJ.SetPoint(i-1, log_bin_center, nevent / dE / exposure);
    gJ.SetPointError(i-1,0,0, n_e1 / dE / exposure, n_e2 / dE / exposure)    



flux = ROOT.TCanvas("flux","flux",900,600)
flux.SetLogy()

gJ.SetMarkerStyle(21);
gJ.SetTitle(";log_{10}(E/eV);J(E) (eV^{-1}m^{-2}sr^{-1}s^{-1})");
gJ.GetYaxis().SetRangeUser(1e-34,2e-30)
gJ.GetXaxis().SetRangeUser(18,20.5)
gJ.Draw("a,e1p");


spec = ROOT.TCanvas("spec","spec",900,600)
spec.SetLogy()

gE3J.SetMarkerStyle(21);
gE3J.SetTitle(";log_{10}(E/eV);E^{3}*J(E) / 1x10^{24} (eV^{2}m^{-2}sr^{-1}s^{-1})");
gE3J.GetYaxis().SetRangeUser(0.1,10.0)
gE3J.GetXaxis().SetRangeUser(18,20.5)
gE3J.Draw("a,e1p")

# Adding the number of events over the points
t = ROOT.TText()
t.SetTextSize(0.02)
t.SetTextFont(42)
for i in range(1,nxbins_md):
    log_bin_center = hEdata_md.GetBinCenter(i)
    log_bin_width = hEdata_md.GetBinWidth(i)
    log_bin_lower_edge = (log_bin_center - log_bin_width / 2.0)
    log_bin_upper_edge = (log_bin_center + log_bin_width / 2.0)
    bin_center = pow(10,log_bin_center)
    bin_lower_edge = pow(10,log_bin_lower_edge)
    bin_upper_edge = pow(10,log_bin_upper_edge)
    Ec3 = pow(bin_center,3)/1e24
    dE = (bin_upper_edge - bin_lower_edge)
    exposure_md = (hExposure_md.GetBinContent(i+10))
    exposure_br = (hExposure_br.GetBinContent(i+10))
    exposure    = exposure_md + exposure_br                 
    nevent_md = hEdata_md.GetBinContent(i)
    nevent_br = hEdata_br.GetBinContent(i)
    nevent    = nevent_md + nevent_br
    if nevent < 20.0 and nevent > 0.0:
        n_e1 = nevent - fldc_err.CalculateLowerLimit(nevent,0.0)
        n_e2 = fldc_err.CalculateUpperLimit(nevent,0.0) - nevent
    elif nevent >= 20:
        n_e1 = np.sqrt(nevent)
        n_e2 = n_e1
    else:
        n_e1 = 0
        n_e2 = 0
    if nevent > 0:
        t.DrawText(log_bin_lower_edge,Ec3 * (nevent + n_e2) / dE / exposure + 0.1,"{0}".format(nevent))
    

# Geting the TA combined spectrum data from root file
ta_comp = ROOT.TFile("{0}/ta_combined_icrc2019.root".format(os.environ['TAX4PROC']),"r")
compspec = ta_comp.Get("compspec")

# Setting variables for TGraph
nentries = compspec.GetEntries()
BC = []
E3J = []
kRed = 2
    
# Fill arrays of values for a TGraph
for i in range(0,nentries):
    compspec.GetEntry(i)
    BC.append(compspec.BinCenter)
    E = pow(10,compspec.BinCenter)
    J = compspec.J
    E3 = pow(E,3)/1e24
    E3J.append(E3*J)

# Setting up TGraph
icrc_compspec = ROOT.TGraph(nentries,np.array(BC),np.array(E3J))
icrc_compspec.SetMarkerColor(kRed);
icrc_compspec.SetMarkerStyle(20);
icrc_compspec.Draw("P");
ROOT.gPad.Update()


# Prompt exit
print('\nPress enter to continue.\n')
os.system('read gonext')
