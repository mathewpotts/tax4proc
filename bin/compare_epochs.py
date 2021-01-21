#!/usr/bin/env python3

# Import Libs
import ROOT
import TDSTio
import TLPY
import os
import argparse

# Call in external C script with quality cut functions
ROOT.gROOT.LoadMacro("hybrid/qualct.C")

# PASS5 FILE PATHS FOR EACH EPOCH
epoch_path = os.environ['MDTAX4_DATA_PASS5'] + '/processing/root_all_files/'
epoch0_mc  = '/home/potts/data/mdtax4/mc/epoch0/epoch0.ps5.drmdpcgf.root'
epoch0     = epoch_path + 'epoch0_20180608-20190409.mdtax4ps5.ps2g.drmdpcgf.root'
epoch1     = epoch_path + 'epoch1_20190424-20190607.mdtax4ps5.ps2g.drmdpcgf.root'
epoch1_mc  = '/home/potts/data/mdtax4/mc/epoch1/epoch1.ps5.drmdpcgf.root'
epoch1_2   = epoch_path + 'epoch0-1.mdtax4.ps5.root'
epoch2     = epoch_path + '20190625-20201206.mdtax4.ps5.ps2g.drmdpcgf.root'
epoch2_mc  = '/home/potts/data/mdtax4/mc/processing/root_all_files/pr.1e17-1e21.20190625-20201206.mdtax4.ps5.drmdpcgf.root'


############# 40 bins e0 vs e0mc Data-MC comparison n_{z} ###################
file = ROOT.TFile(epoch0,'r')
taTree_e0 = file.Get('taTree')
hepoch0 = ROOT.TH1D('hepoch0','Epoch 0 n_{z} Data-MC comparison;n_{z};',40,-1,1)
taTree_e0.Draw('stpln.n_ampwt[6][2]>>hepoch0','etrack.qualct==1','goff')

file1 = ROOT.TFile(epoch0_mc,'r')
taTree_e0mc = file1.Get('taTree')
hepoch0mc = ROOT.TH1D('hepoch0mc','hepoch0mc',40,-1,1)
taTree_e0mc.Draw('stpln.n_ampwt[6][2]>>hepoch0mc','weight(log10(mc04->energy))*(etrack.qualct==1)','goff')
norm = hepoch0.Integral()/hepoch0mc.Integral()
hepoch0mc.Scale(norm)

c2e0rp = ROOT.TCanvas("c2e0rp","c2e0rp",800,600)
rpe0 = ROOT.TRatioPlot(hepoch0,hepoch0mc,"divsym")
rpe0.SetH1DrawOpt("E0")
rpe0.SetH2DrawOpt("HIST")
c2e0rp.cd()
rpe0.Draw()
rpe0.GetLowerRefYaxis().SetRangeUser(0.1,3)
rpe0.GetLowerRefYaxis().SetLabelSize(0.025)
rpe0.GetLowerRefGraph().Fit("pol1","W")

ROOT.gPad.BuildLegend()
ROOT.gPad.Update()

# Prompt exit
print('\nPress enter to continue.\n')
os.system('read gonext')

############# 40 bins e1 vs e1mc Data-MC comparison n_{z} ###################
file = ROOT.TFile(epoch1,'r')
taTree_e1 = file.Get('taTree')
hepoch1 = ROOT.TH1D('hepoch1','Epoch 1 n_{z} Data-MC comparison;n_{z};',40,-1,1)
taTree_e1.Draw('stpln.n_ampwt[6][2]>>hepoch1','etrack.qualct==1')

file1 = ROOT.TFile(epoch1_mc,'r')
taTree_e1mc = file1.Get('taTree')
hepoch1mc = ROOT.TH1D('hepoch1mc','hepoch1mc',40,-1,1)
taTree_e1mc.Draw('stpln.n_ampwt[6][2]>>hepoch1mc','weight(log10(mc04->energy))*(etrack.qualct==1)','p sames')
norm = hepoch1.Integral()/hepoch1mc.Integral()
hepoch1mc.Scale(norm)

c2e1rp = ROOT.TCanvas("c2e1rp","c2e1rp",800,600)
rpe1 = ROOT.TRatioPlot(hepoch1,hepoch1mc,"divsym")
rpe1.SetH1DrawOpt("E0")
rpe1.SetH2DrawOpt("HIST")
c2e1rp.cd()
rpe1.Draw()
rpe1.GetLowerRefYaxis().SetRangeUser(0.1,3)
rpe1.GetLowerRefYaxis().SetLabelSize(0.025)
rpe1.GetLowerRefGraph().Fit("pol1","W")

ROOT.gPad.BuildLegend()
ROOT.gPad.Update()

# Prompt exit
print('\nPress enter to continue.\n')
os.system('read gonext')

############# 40 bins e2 vs e2mc Data-MC comparison n_{z} ###################
file = ROOT.TFile(epoch2,'r')
taTree_e2 = file.Get('taTree')
hepoch2 = ROOT.TH1D('hepoch2','Epoch 2 n_{z} Data-MC comparison;n_{z};',40,-1,1)
taTree_e2.Draw('stpln.n_ampwt[6][2]>>hepoch2','etrack.qualct==1')

file1 = ROOT.TFile(epoch2_mc,'r')
taTree_e2mc = file1.Get('taTree')
hepoch2mc = ROOT.TH1D('hepoch2mc','hepoch2mc',40,-1,1)
taTree_e2mc.Draw('stpln.n_ampwt[6][2]>>hepoch2mc','weight(log10(mc04->energy))*(etrack.qualct==1)','p sames')
norm = hepoch2.Integral()/hepoch2mc.Integral()
hepoch2mc.Scale(norm)

c2e2rp = ROOT.TCanvas("c2e2rp","c2e2rp",800,600)
rpe2 = ROOT.TRatioPlot(hepoch2,hepoch2mc,"divsym")
rpe2.SetH1DrawOpt("E0")
rpe2.SetH2DrawOpt("HIST")
c2e2rp.cd()
rpe2.Draw()
rpe2.GetLowerRefYaxis().SetRangeUser(0.1,3)
rpe2.GetLowerRefYaxis().SetLabelSize(0.025)
rpe2.GetLowerRefGraph().Fit("pol1","W")

ROOT.gPad.BuildLegend()
ROOT.gPad.Update()

# Prompt exit
print('\nPress enter to continue.\n')
os.system('read gonext')

############ 30 bins e0 vs e0mc Data-MC comparison cos^{-1}(n_{z}) ############
file = ROOT.TFile(epoch0,'r')
taTree_e0 = file.Get('taTree')
hepoch0 = ROOT.TH1D('hepoch0','Epoch 0 cos^{-1}(n_{z}) Data-MC comparison;cos^{-1}(n_{z}) degrees;',30,-90,90)
taTree_e0.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2]))- 90>>hepoch0','etrack.qualct==1','goff')

file1 = ROOT.TFile(epoch0_mc,'r')
taTree_e0mc = file1.Get('taTree')
hepoch0mc = ROOT.TH1D('hepoch0mc','hepoch0mc',30,-90,90)
taTree_e0mc.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2]))- 90>>hepoch0mc','weight(log10(mc04->energy))*(etrack.qualct==1)','goff')
norm = hepoch0.Integral()/hepoch0mc.Integral()
hepoch0mc.Scale(norm)

c2e0rp = ROOT.TCanvas("c2e0rp","c2e0rp",800,600)
rpe0 = ROOT.TRatioPlot(hepoch0,hepoch0mc,"divsym")
rpe0.SetH1DrawOpt("E0")
rpe0.SetH2DrawOpt("HIST")
c2e0rp.cd()
rpe0.Draw()
rpe0.GetLowerRefYaxis().SetRangeUser(0.1,3)
rpe0.GetLowerRefYaxis().SetLabelSize(0.025)
rpe0.GetLowerRefGraph().Fit("pol1","W")

ROOT.gPad.BuildLegend()
ROOT.gPad.Update()

# Prompt exit
print('\nPress enter to continue.\n')
os.system('read gonext')
 
############ 30 bins e1 vs e1mc Data-MC comparison cos^{-1}(n_{z}) ############
file = ROOT.TFile(epoch1,'r')
taTree_e1 = file.Get('taTree')
hepoch1 = ROOT.TH1D('hepoch1','Epoch 1 cos^{-1}(n_{z}) Data-MC comparison;cos^{-1}(n_{z}) degrees;',30,-90,90)
taTree_e1.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2]))- 90>>hepoch1','etrack.qualct==1')

file1 = ROOT.TFile(epoch1_mc,'r')
taTree_e1mc = file1.Get('taTree')
hepoch1mc = ROOT.TH1D('hepoch1mc','hepoch1mc',30,-90,90)
taTree_e1mc.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2]))- 90>>hepoch1mc','weight(log10(mc04->energy))*(etrack.qualct==1)','p sames')
norm = hepoch1.Integral()/hepoch1mc.Integral()
hepoch1mc.Scale(norm)

c2e1rp = ROOT.TCanvas("c2e1rp","c2e1rp",800,600)
rpe1 = ROOT.TRatioPlot(hepoch1,hepoch1mc,"divsym")
rpe1.SetH1DrawOpt("E0")
rpe1.SetH2DrawOpt("HIST")
c2e1rp.cd()
rpe1.Draw()
rpe1.GetLowerRefYaxis().SetRangeUser(0.1,3)
rpe1.GetLowerRefYaxis().SetLabelSize(0.025)
rpe1.GetLowerRefGraph().Fit("pol1","W")

ROOT.gPad.BuildLegend()
ROOT.gPad.Update()

# Prompt exit
print('\nPress enter to continue.\n')
os.system('read gonext')

############ 30 bins e2 vs e2mc Data-MC comparison cos^{-1}(n_{z}) ############
file = ROOT.TFile(epoch2,'r')
taTree_e2 = file.Get('taTree')
hepoch2 = ROOT.TH1D('hepoch2','Epoch 2 cos^{-1}(n_{z}) Data-MC comparison;cos^{-1}(n_{z}) degrees;',30,-90,90)
taTree_e2.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2]))- 90>>hepoch2','etrack.qualct==1')

file1 = ROOT.TFile(epoch2_mc,'r')
taTree_e2mc = file1.Get('taTree')
hepoch2mc = ROOT.TH1D('hepoch2mc','hepoch2mc',30,-90,90)
taTree_e2mc.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2]))- 90>>hepoch2mc','weight(log10(mc04->energy))*(etrack.qualct==1)','p sames')
norm = hepoch2.Integral()/hepoch2mc.Integral()
hepoch2mc.Scale(norm)

c2e2rp = ROOT.TCanvas("c2e2rp","c2e2rp",800,600)
rpe2 = ROOT.TRatioPlot(hepoch2,hepoch2mc,"divsym")
rpe2.SetH1DrawOpt("E0")
rpe2.SetH2DrawOpt("HIST")
c2e2rp.cd()
rpe2.Draw()
rpe2.GetLowerRefYaxis().SetRangeUser(0.1,3)
rpe2.GetLowerRefYaxis().SetLabelSize(0.025)
rpe2.GetLowerRefGraph().Fit("pol1","W")

ROOT.gPad.BuildLegend()
ROOT.gPad.Update()

# Prompt exit
print('\nPress enter to continue.\n')
os.system('read gonext')

############ 30 bins all epochs ############
file = ROOT.TFile(epoch0_mc,'r')
taTree_e0 = file.Get('taTree')
nentries_e0 = taTree_e0.GetEntries()
hepoch0 = ROOT.TH1D('hepoch0mc','All epochs MC cos^{-1}(n_{z}) comparison;cos^{-1}(n_{z}) degrees;',30,-90,90)
taTree_e0.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2]))- 90>>hepoch0mc','weight(log10(mc04->energy))*(etrack.qualct==1)')

file1 = ROOT.TFile(epoch1_mc,'r')
taTree_e1 = file1.Get('taTree')
hepoch1 = ROOT.TH1D('hepoch1mc','hepoch1mc',30,-90,90)
taTree_e1.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2]))- 90>>hepoch1mc','weight(log10(mc04->energy))*(etrack.qualct==1)','p sames')
norm = hepoch0.Integral()/hepoch1.Integral()
hepoch1.Scale(norm)

file2 = ROOT.TFile(epoch2_mc,'r')
taTree_e2 = file2.Get('taTree')
hepoch2 = ROOT.TH1D('hepoch2mc','hepoch2mc',30,-90,90)
taTree_e2.Draw('(180/TMath::Pi())*(TMath::ACos(stpln.n_ampwt[6][2]))- 90>>hepoch2mc','weight(log10(mc04->energy))*(etrack.qualct==1)','p sames')
norm = hepoch0.Integral()/hepoch2.Integral()
hepoch2.Scale(norm)
ROOT.gPad.BuildLegend()
ROOT.gPad.Update()

# Prompt exit
print('\nPress enter to continue.\n')
os.system('read gonext')
